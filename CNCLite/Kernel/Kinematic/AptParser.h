#ifndef APTPARSER_H
#define APTPARSER_H
#include <QFile>
#include <QString>
#include <QRegExp>
#include "../../ThirdParty/Eigen/Dense"

namespace CNCLite {
class AptParser
{
public:
    AptParser(const QString& fn);
    ~AptParser();
public:
    Eigen::VectorXd readLine();
    void readAll(Eigen::MatrixXd& data, uint interval=1);

private:
    QFile aptFile;
    bool isFileValid;
    // The regular expression for valid position and orientation data in APT files.
    QRegExp reFormat;
};
} /// End of namespace CNCLite.

#endif // APTPARSER_H
