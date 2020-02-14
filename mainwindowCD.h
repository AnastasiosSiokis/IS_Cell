#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "IS_cell.h"

#include "qcustomplotCD.h"

QMap<double, QColor> AdaptColorScale(bool isCheckedTCR, bool isCheckedLFA1, bool isCheckedpMHC, bool isCheckedICAM1, bool isCheckedTM, bool isCheckedLI, bool isCheckedCD45,
                                     bool isCheckedCD2, bool isCheckedCD48, bool isCheckedCD2_CD48, bool isCheckedCD28, bool isCheckedCD80, bool isCheckedCD28_CD80,
                                     bool isCheckedPD1, bool isCheckedPDL1, bool isCheckedPD1_PDL1);

class Plotting : public QWidget{
        Q_OBJECT
public:
    QCustomPlot* customPlot;
    Plotting(QWidget *parent, std::string title = std::string(""));
    QCPColorMap *colorMap;
    QCPColorScale *colorScale;
    QCPColorGradient *myGrad;
    void PlotLattice(lattice *Cell);
    void PlotLattice(doubleLattice *Cell);
    void ChangeColors(QMap<double, QColor> colStops);
};


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    Plotting* plotForTcell;
    Plotting* plotForAPC;
    Plotting* plotForColocalization;
    Plotting* plotForTracking;
    Plotting* plotForActin;

    #ifdef USE_DENS_LATTICES
    Plotting* plotForDensT;
    Plotting* plotForDensM;
    Plotting* plotForDensTM;
    Plotting* plotForCellKD;
    #endif


    void PlotLatt(lattice* Cell);
    void RunSimulation();

    simulation* gamoto;
    bool isRunning;
    double t;

public slots:
    void restart();
    void runOrStop();
    void TcellToFig(QString fileNameTcell = QString(""));
    void APCToFig(QString fileNameAPC = QString(""));
    void ColocalizationToFig(QString fileNameActin = QString(""));
    void TrackingToFig(QString fileNameTracking = QString(""));
    void ActinToFig(QString fileNameTracking = QString(""));
    void DensTToFig(QString fileNameDensT = QString(""));
    void DensMToFig(QString fileNameDensM = QString(""));
    void DensTMToFig(QString fileNameDensTM = QString(""));
    void CellKDToFig(QString fileNameCellKD = QString(""));
    void MakeAllFig();

    void TCRcheck(int);
    void LFA1check(int);
    void pMHCcheck(int);
    void ICAM1check(int);
    void TMcheck(int);
    void LIcheck(int);
    void CD2_CD48check(int);
    void CD28_CD80check(int);
    void PD1_PDL1check(int);
    void CD45check(int);
    void CD2check(int);
    void CD48check(int);
    void CD28check(int);
    void CD80check(int);
    void PD1check(int);
    void PDL1check(int);
    void Koff_LI(double val1);
    void Koff_TM(double val2);
    void Kon_LI(double val3);
    void Kon_TM(double val4);

    void updateColors();

    void getDefaultValues();

    void getCursorCoord(QMouseEvent *event);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H

