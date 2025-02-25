#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class QMenu;
class QAction;
class CfdThread;
class QProcess;
class QSplitter;
class Terminal;
class Panel;

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    void closeEvent( QCloseEvent * ev ) override;
protected:
    bool eventFilter(QObject *obj, QEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;
private:
    void triggerNew();
    void runCFD();
    void runMPI();
    void runTerminal();
private:
    Ui::MainWindow *ui;
    QMenu * menuFile;
    QAction * actNew;
    QAction * actRun;
    QAction * actRunMPI;
    QAction * actTerminal;
private:
    CfdThread *cfdThread;
    Terminal * terminal = nullptr;
    Panel * panel = nullptr;
};
#endif // MAINWINDOW_H
