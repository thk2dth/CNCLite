#include "Trajectory.h"

using namespace Eigen;
using namespace CNCLite;

Trajectory::Trajectory()
{
    isPlanned = false;
    duration = 0.0;
    trajType = UNDEFINED;
    law.resize(maxNumLaw);
    sa = 0.0;
    sc = 0.0;
    sd = 0.0;
    ta = 0.0;
    tc = 0.0;
    td = 0.0;
}

Trajectory::Trajectory(double vs, double ve, const VectorXd &kineConst, Path *path,
                       TimeLaw::LawType type, double us, double ue)
{
    initialize(vs, ve, kineConst, path, type, us, ue);
}

Trajectory::~Trajectory()
{
    for(unsigned int i=0; i<maxNumLaw; i++)
    {
        if (law.at(i) != nullptr)
        {
            delete law[i];
            law[i] = nullptr;
        }
    }
}

void Trajectory::initialize(double vs, double ve, const VectorXd &kineConst, Path *path,
                            TimeLaw::LawType type, double us, double ue)
{
    this->vs = vs;
    this->ve = ve;
    kinematicConstraint = kineConst;
    this->path = path;
    this->us = us;
    this->ue = ue;
    uc = us;
    dc = 0.0;
    isPlanned = false;
    duration = 0.0;
    trajType = UNDEFINED;
    law.resize(maxNumLaw);
    if (fabs(us-path->getBeginParameter())<=EPS_NUM &&
            fabs(ue-path->getEndParameter())<=EPS_NUM )
        length = path->getLength();
    else
        length = path->calculateLengthInterval(us, ue);
    sa = 0.0;
    sc = 0.0;
    sd = 0.0;
    ta = 0.0;
    tc = 0.0;
    td = 0.0;
    for (uint8_t i=0; i<maxNumLaw; i++)
        law[i] = nullptr;
    lawType = type;
}

bool Trajectory::plan()
{
    unsigned int order = kinematicConstraint.size(); // order of the time law.
    if (law[0] == nullptr)
    {
        /// We do not need to reallocate memory for law.
        switch (lawType)
        {
        case TimeLaw::POLYNOMIAL:
            for (int8_t i=0; i<maxNumLaw; i++)
            {
                switch (order)
                {
                case 2:
                    law[i] = new AccelerationBounded();
                    break;
                case 3:
                    law[i] = new JerkBounded();
                    break;
                case 4:
                    law[i] = new SnapBounded();
                    break;
                default:
                    law[i] = new SnapBounded();
                }
            }
            break;
        case TimeLaw::SINE_SERIES:
            for (int8_t i=0; i<maxNumLaw; i++)
                law[i] = new JerkSineSeries();
            break;
        case TimeLaw::SINE:
            for (int8_t i=0; i<maxNumLaw; i++)
                law[i] = new JerkSine();
            break;
        }
    }
    law[0]->initialize(kinematicConstraint, true);
    law[1]->initialize(kinematicConstraint, true);
    law[2]->initialize(kinematicConstraint, false);
    double s0; // minimum distance required to accelerate from vs to ve.
    bool isEndGreatBegin = (vs <= ve);
    if (vs == ve)
        s0 = 0.0;
    else if (vs < ve)
        s0 = law[0]->lawVV(vs, ve);
    else
        s0 = law[2]->lawVV(vs, ve);
    if (length <= s0)
    {
        /// Only one time law is needed.
        if (isEndGreatBegin)
        {
            ve = law[0]->lawVD(vs, length, true);
            vf = ve;
            trajType = ACC;
            sa = length;
            ta = law[0]->getDuration();
            duration = law[0]->getDuration();
        }
        else
        {
            vs = law[2]->lawVD(ve, length, false);
            vf = vs;
            trajType = DEC;
            sd = length;
            td = law[2]->getDuration();
            duration = law[2]->getDuration();
        }
        isPlanned = true;
        return true;
    }
    else
    {
        double vm = kinematicConstraint(0); /// alias of maximum velocity.
        assert(vm > 0);
        double high = vm;
        sa = law[0]->lawVV(vs, vm);
        sd = law[2]->lawVV(vm, ve);
        if (length < sa + sd)
        {
            /// Binary search.
            double low = std::fmax(vs, ve); // search interval [low, high].
            /// In case fmax(vs,ve)=vm, do lawVV() first.
            do
            {
                vf = 0.5 * (high+low);
                sa = law[0]->lawVV(vs, vf);
                sd = law[2]->lawVV(vf, ve); // DEC phase.
                if(sa+sd <= length)
                {
                    low = vf;
                }
                else
                {
                    high = vf;
                }
            }while (fabs(sa + sd - length) > EPS_CAD);
            isPlanned = true;
            trajType = ACC_DEC;
            ta = law[0]->getDuration();
            td = law[2]->getDuration();
            duration = ta + td;
            law[0]->setIsAcc(true);
            law[2]->setIsAcc(false);
            sc = 0.0;
            tc = 0.0;
            return false;
        }
        else
        {
            isPlanned = true;
            trajType = ACC_CRU_DEC;
            vf = vm;
            sc = length - sa - sd;
            law[1]->lawCruise(vf, sc);
            ta = law[0]->getDuration();
            td = law[2]->getDuration();
            tc = law[1]->getDuration();
            duration = ta + td + tc;
            return false;
        }
    }
    return false;
}

bool Trajectory::synchronize(double tg)
{
    if (!isPlanned)
        plan();
    if (tg <= duration)
    {
        std::cout << "Cannot synchronize the traversing time to a shorter duration.\n";
        return true;
    }
    double T1, T2, T3;
    if (ve >= vs)
    {
        /// Eq. (3).
        /// If ve=0, T1 will be undefined.
        if (ve == 0)
            T1 = std::numeric_limits<double>::max();
        else
            T1 = (length + 0.5 * (ve -vs)*ta ) / ve;
        if (tg <= T1)
        {
            /// Eq. (4).
            vf = (tc+0.5*(ta+td)) / (tg-0.5*(ta+td)) * vf;
            law[0]->stretchByV(vf, true);
            law[2]->stretchByV(vf, false);
            sa = 0.5 * (vf+vs) * ta;
            sd = 0.5 * (vf+ve) * td;
            sc = length - sa - sd;
            law[1]->lawCruise(vf, sc);
            tc = law[1]->getDuration(); // tc is increased.
            switch(trajType)
            {
            case ACC:
                /// This situation normally does not occur.
                trajType = ACC_CRU; break;
            case ACC_DEC:
            case ACC_CRU_DEC:
                trajType = ACC_CRU_DEC; break;
            default:
                std::cerr<<"Uncaptured situations during time synchronization.\n"; break;
            }
            duration = ta + tc + td;
            /// assert(std::fabs(tg-duration) < EPS_NUM);
            return false;
        }
        else
        {
            /// Eq. (6).
            T2 = 2*length / (vs+ve);
            if (tg <= T2)
            {
                ta = 2*(ve*tg-length) / (ve-vs); // ta is extended to Tl.
                law[0]->stretchByV(ve, true);
                law[0]->stretchByTime(ta);
                sa = law[0]->getDispm();
                sc = length - sa;
                sd = 0;
                td = 0.0;
                law[1]->lawCruise(ve, sc);
                tc = law[1]->getDuration(); // tc is increased.
                law[2]->reset(); // No DEC phase.
                trajType = ACC_CRU;
                duration = tg;
                return false;
            }
            else
            {
                if (vs == 0)
                    T3 = std::numeric_limits<double>::max();
                else
                    T3 = length / vs;
                if (tg <= T3)
                {
                    vf = ve*length / (tg*(vs+ve)-length);
                    law[0]->stretchByV(vf, true);
                    law[0]->stretchByTime(T2);
                    sa = law[0]->getDispm();
                    ta = T2;
                    sc = length - sa;
                    sd = 0;
                    td = 0.0;
                    law[1]->lawCruise(vf, sc);
                    tc = law[1]->getDuration();
                    law[2]->reset(); // No DEC phase.
                    trajType = ACC_CRU;
                    duration = tg;
                    ve = vf;
                    return true;
                }
                else
                {
                    /// Before stretching, the third time law should be re-planned.
                    law[0]->lawVV(vs, 0.0); // Decelerate from vs to 0.
                    td = law[0]->getDuration();

                    /// Eq. (12)
                    vf = (2*length-vs*td) / (2*tg-td);
                    assert(vf >= 0);
                    law[0]->stretchByV(vf, true); // Increase the end velocity.
                    sd = law[0]->getDispm();
                    sc = length - sd;
                    sa = 0.0;
                    ta = 0.0;
                    law[1]->lawCruise(vf, sc);
                    tc = law[1]->getDuration();
                    law[2]->reset();
                    trajType = DEC_CRU;
                    ve = vf;
                    duration = tg;
                    return true;
                }
            }
        }
    }
    else /// vs > ve
    {
        /// Eq. (3).
        if (vs == 0)
            T1 = std::numeric_limits<double>::max();
        else
            T1 = (length + 0.5 * (vs -ve)*td ) / vs;
        if (tg <= T1)
        {
            /// Eq. (4).
            vf = (tc+0.5*(ta+td)) / (tc+0.5*(ta+td)+tg-duration) * vf;
            law[0]->stretchByV(vf, true);
            law[2]->stretchByV(vf, false);
            sa = 0.5*(vf + vs)*ta;
            sd = 0.5*(vf+ve)*td;
            sc = length - sa - sd;
            law[1]->lawCruise(vf, sc);
            tc = law[1]->getDuration(); // tc is increased.
            switch(trajType)
            {
            case DEC:
                /// This situation normally does not occur.
                trajType = CRU_DEC; break;
            case ACC_DEC:
            case ACC_CRU_DEC:
                trajType = ACC_CRU_DEC; break;
            default:
                std::cerr<<"Uncaptured situations during time synchronization.\n"; break;
            }
            duration = tg;
            return false;
        }
        else
        {
            /// Eq. (6).
            T2 = 2*length / (vs+ve);
            if (tg <= T2)
            {
                td = 2*(vs*tg-length) / (vs-ve);
                law[2]->stretchByV(vs, false);
                law[2]->stretchByTime(td);
                sd = law[2]->getDispm();
                sc = length - sd;
                sa = 0;
                ta = 0;
                law[1]->lawCruise(vs, sc);
                tc = law[1]->getDuration();
                law[0]->reset(); // No ACC phase.
                trajType = CRU_DEC;
                duration = tg;
                return false;
            }
            else
            {
                if (ve == 0)
                    T3 = std::numeric_limits<double>::max();
                else
                    T3 = length / ve;
                if (tg <= T3)
                {
                    vf = vs*length / (tg*(vs+ve)-length);
                    law[2]->stretchByV(vf, false);
                    law[2]->stretchByTime(T2);
                    sd = law[2]->getDispm();
                    td = T2;
                    sc = length - sd;
                    sa = 0;
                    ta = 0;
                    law[1]->lawCruise(vf, sc);
                    tc = law[1]->getDuration();
                    law[0]->reset(); // No DEC phase.
                    trajType = CRU_DEC;
                    duration = tg;
                    vs = vf;
                    return true;
                }
                else
                {
                    /// Before stretching, the third time law should be re-planned.
                    law[2]->lawVV(0.0, ve);
                    ta = law[2]->getDuration();

                    /// Eq. (12)
                    vf = (2*length-ve*ta) / (2*tg-ta);
                    law[2]->stretchByV(vf, false); // Increase the start velocity.
                    sa = law[2]->getDispm();
                    sc = length - sa;
                    sd = 0.0;
                    td = 0;
                    law[1]->lawCruise(vf, sc);
                    tc = law[1]->getDuration();
                    law[0]->reset();
                    trajType = CRU_ACC;
                    vs = vf;
                    duration = tg;
                    return true;
                }
            }
        }
    }
}


double Trajectory::displacement(double t) const
{
    double t_bgn = 0.0;
    double t_end = 0.0;
    double accum = 0.0;
    if (t >= duration)
        return length;
    else if (t <= 0)
        return 0.0;
    else
    {
        for (unsigned int i=0; i<maxNumLaw; i++)
        {
            TimeLaw* tl = law.at(i);
            t_bgn = t_end;
            t_end += tl->getDuration();
            if (t<=t_end)
            {
                return tl->displacement(t-t_bgn) + accum;
                break;
            }
            else
            {
                accum += tl->getDispm(); // accumulative displacement.
                continue;
            }
        }
    }
}

double Trajectory::velocity(double t) const
{
    double t_bgn = 0.0;
    double t_end = 0.0;
    if (t >= duration)
        return ve;
    else if (t <= 0)
        return vs;
    else
    {
        for (unsigned int i=0; i<maxNumLaw; i++)
        {
            TimeLaw *tl = law.at(i);
            t_bgn = t_end;
            t_end += tl->getDuration();
            if (t<=t_end)
            {
                return tl->velocity(t-t_bgn);
                break;
            }
            else
                continue;
        }
    }
}

double Trajectory::acceleration(double t) const
{
    double t_bgn = 0.0;
    double t_end = 0.0;
    if (t >= duration)
        return 0.0;
    else if (t <= 0)
        return 0.0;
    else
    {
        for (unsigned int i=0; i<maxNumLaw; i++)
        {
            TimeLaw *tl = law.at(i);
            t_bgn = t_end;
            t_end += tl->getDuration();
            if (t<=t_end)
            {
                return tl->acceleration(t-t_bgn);
                break;
            }
            else
                continue;
        }
    }
}

double Trajectory::jerk(double t) const
{
    double t_bgn = 0.0;
    double t_end = 0.0;
    if (t >= duration)
        return 0.0;
    else if (t <= 0)
        return 0.0;
    else
    {
        for (unsigned int i=0; i<maxNumLaw; i++)
        {
            TimeLaw *tl = law.at(i);
            t_bgn = t_end;
            t_end += tl->getDuration();
            if (t<=t_end)
            {
                return tl->jerk(t-t_bgn);
                break;
            }
            else
                continue;
        }
    }
}

double Trajectory::snap(double t) const
{
    double t_bgn = 0.0;
    double t_end = 0.0;
    if (t >= duration)
        return 0.0;
    else if (t <= 0)
        return 0.0;
    else
    {
        for (unsigned int i=0; i<maxNumLaw; i++)
        {
            TimeLaw *tl = law.at(i);
            t_bgn = t_end;
            t_end += tl->getDuration();
            if (t<=t_end)
            {
                return tl->snap(t-t_bgn);
                break;
            }
            else
                continue;
        }
    }
    return 0.0;
}

/// Function axialPosition() will update current parameter uc.
VectorXd Trajectory::axialPosition(double t)
{
    double d = displacement(t);
    double ds = d - dc;
    dc = d;
    uc = path->calculateNextParameter(ds, uc, Path::SECONDRUKUCOM);
    return path->calculatePoint(uc);
}

VectorXd Trajectory::axialVelocity(double t)
{
    double d = displacement(t);
    /// If t remains the same, ds will be 0. uc will not be updated.
    double ds = d - dc;
    dc = d;
    uc = path->calculateNextParameter(ds, uc, Path::SECONDRUKUCOM);
    VectorXd vec = path->calculateDer1(uc).normalized();
    double f = velocity(t);
    return vec * f;
}

VectorXd Trajectory::axialPosition(double ds, double u)
{
    double un = path->calculateNextParameter(ds, u, Path::SECONDRUKUCOM);
    return path->calculatePoint(un);
}

VectorXd Trajectory::axialVelocity(double ds, double u, double f)
{
    double un = path->calculateNextParameter(ds, u, Path::SECONDRUKUCOM);
    VectorXd vec = path->calculateDer1(un).normalized();
    return vec * f;
}

VectorXd Trajectory::getKinematicConstraint() const
{
    return kinematicConstraint;
}

double Trajectory::getDuration() const
{
    return duration;
}

Path* Trajectory::getPath() const
{
    return path;
}

double Trajectory::getVs() const
{
    return vs;
}

double Trajectory::getVe() const
{
    return ve;
}

double Trajectory::getVf() const
{
    return vf;
}

double Trajectory::getLength() const
{
    return length;
}

void Trajectory::setVs(double v)
{
    vs = v;
}

void Trajectory::setVe(double v)
{
    ve = v;
}

double Trajectory::getUc() const
{
    return uc;
}

Trajectory::TrajType Trajectory::getTrajType() const
{
    return trajType;
}

double Trajectory::getUs() const
{
    return us;
}

double Trajectory::getUe() const
{
    return ue;
}
