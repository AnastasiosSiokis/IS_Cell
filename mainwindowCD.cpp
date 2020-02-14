#include "mainwindowCD.h"
#include "ui_mainwindowCD.h"
#include <string>
#include <cmath> //(phi)
#include <QApplication>
#include <iostream>
#include <QMap>
#include <QFileDialog>
//#include <QtSvg>

QMap<double, QColor> AdaptColorScale(bool isCheckedTCR, bool isCheckedLFA1, bool isCheckedpMHC, bool isCheckedICAM1, bool isCheckedTM, bool isCheckedLI, bool isCheckedCD45,
                                     bool isCheckedCD2, bool isCheckedCD48, bool isCheckedCD2_CD48, bool isCheckedCD28, bool isCheckedCD80, bool isCheckedCD28_CD80,
                                     bool isCheckedPD1, bool isCheckedPDL1, bool isCheckedPD1_PDL1){
   // myGrad->setLevelCount(21);                  // number of color levels in the gradient
    QMap<double, QColor> colStops;
#define coef 10./27.0
    colStops.insert(0.00 * coef,  QColor(Qt::black)); // nothing
    colStops.insert(0.05 * coef,  QColor(Qt::black));// border
    colStops.insert(0.10 * coef,  QColor(Qt::black));
    colStops.insert(0.15 * coef,  QColor((isCheckedTCR)? Qt::yellow : Qt::black)); // TCR
    colStops.insert(0.20 * coef,  QColor((isCheckedTCR)? Qt::yellow : Qt::black));
    colStops.insert(0.25 * coef,  QColor((isCheckedLFA1)? Qt::cyan : Qt::black)); // LFA1 is Orange --> 255,128,0 QColor(255,255,255)
    colStops.insert(0.30 * coef,  QColor((isCheckedLFA1)? Qt::cyan : Qt::black));
    colStops.insert(0.35 * coef,  QColor((isCheckedpMHC)? Qt::cyan : Qt::black)); // pMHC
    colStops.insert(0.40 * coef,  QColor((isCheckedpMHC)? Qt::cyan : Qt::black));
    colStops.insert(0.45 * coef,  QColor((isCheckedICAM1)? Qt::magenta : Qt::black)); // ICAM1
    colStops.insert(0.50 * coef,  QColor((isCheckedICAM1)? Qt::magenta : Qt::black));
    colStops.insert(0.55 * coef,  QColor((isCheckedTM)? Qt::green : Qt::black)); // TCR-pMHC
    colStops.insert(0.60 * coef,  QColor((isCheckedTM)? Qt::green : Qt::black));
    colStops.insert(0.65 * coef,  QColor((isCheckedLI)? Qt::red : Qt::black)); // LFA1-ICAM1 (220,20,60)
    colStops.insert(0.70 * coef,  QColor((isCheckedLI)? Qt::red : Qt::black));
    colStops.insert(0.75 * coef,  QColor(0,191,255)); // Nucleation Points
    colStops.insert(0.80 * coef,  QColor(0,191,255));
    colStops.insert(0.85 * coef,  QColor(Qt::red)); // Nucleation Polymerized Points
    colStops.insert(0.90 * coef,  QColor(Qt::red));
    colStops.insert(0.95 * coef,  QColor(Qt::yellow)); // Polymerization Points
    colStops.insert(1.00 * coef,  QColor(Qt::yellow));
    colStops.insert(1.05 * coef,  QColor((isCheckedCD45)? Qt::white : Qt::black)); // CD45
    colStops.insert(1.10 * coef,  QColor((isCheckedCD45)? Qt::white : Qt::black));
    colStops.insert(1.15 * coef,  QColor(Qt::black)); // Empty Actin Points
    colStops.insert(1.20 * coef,  QColor(Qt::black));
    colStops.insert(1.25 * coef,  QColor((isCheckedCD2)? Qt::cyan : Qt::black)); // CD2
    colStops.insert(1.30 * coef,  QColor((isCheckedCD2)? Qt::cyan : Qt::black));
    colStops.insert(1.35 * coef,  QColor((isCheckedCD48)? Qt::yellow : Qt::black)); // CD48
    colStops.insert(1.40 * coef,  QColor((isCheckedCD48)? Qt::yellow : Qt::black));
    colStops.insert(1.45 * coef,  QColor((isCheckedCD2_CD48)? Qt::magenta : Qt::black)); // CD2_CD48
    colStops.insert(1.50 * coef,  QColor((isCheckedCD2_CD48)? Qt::magenta: Qt::black));
    colStops.insert(1.55 * coef,  QColor((isCheckedCD28)? Qt::yellow : Qt::black)); // CD28
    colStops.insert(1.60 * coef,  QColor((isCheckedCD28)? Qt::yellow : Qt::black));
    colStops.insert(1.65 * coef,  QColor((isCheckedCD80)? Qt::blue : Qt::black)); // CD80
    colStops.insert(1.70 * coef,  QColor((isCheckedCD80)? Qt::blue : Qt::black));
    colStops.insert(1.75 * coef,  QColor((isCheckedCD28_CD80)? Qt::cyan : Qt::black)); // CD28_CD80 rgb(0,191,255)
    colStops.insert(1.80 * coef,  QColor((isCheckedCD28_CD80)? Qt::cyan : Qt::black));
    colStops.insert(1.85 * coef,  QColor(Qt::white));
    colStops.insert(1.90 * coef,  QColor(Qt::white));
    colStops.insert(1.95 * coef,  QColor(Qt::white));
    colStops.insert(2.00 * coef,  QColor(Qt::white));
    colStops.insert(2.05 * coef,  QColor(Qt::white));
    colStops.insert(2.10 * coef,  QColor(Qt::white));
    colStops.insert(2.15 * coef,  QColor(Qt::white));
    colStops.insert(2.20 * coef,  QColor(Qt::white));
    colStops.insert(2.25 * coef,  QColor(Qt::white));
    colStops.insert(2.30 * coef,  QColor(Qt::white));
    colStops.insert(2.35 * coef,  QColor(Qt::white));
    colStops.insert(2.40 * coef,  QColor(Qt::white));
    colStops.insert(2.45 * coef,  QColor(Qt::white));
    colStops.insert(2.50 * coef,  QColor(Qt::white));
    return colStops;
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow){
    ui->setupUi(this);

    QDir::setCurrent(QCoreApplication::applicationDirPath());
    ui->scrollAreaWidgetContents->resize(1000, 1000);
    ui->scrollArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    plotForTcell = new Plotting(ui->widget, std::string("T Cell"));    // Tcell
    plotForAPC = new Plotting(ui->widget_2, std::string("APC")); // APC
    plotForColocalization = new Plotting(ui->widget_3, std::string("Colocalization")); // Colocalizatio
    plotForTracking = new Plotting(ui->widget_4, std::string("Tracking T cell & APC")); // Actin
    plotForActin = new Plotting(ui->widget_5, std::string("Actin")); // Colocalization

    connect(plotForTcell->customPlot, SIGNAL(mouseRelease(QMouseEvent*)), this,SLOT(getCursorCoord(QMouseEvent*)));
    connect(plotForAPC->customPlot, SIGNAL(mouseRelease(QMouseEvent*)), this,SLOT(getCursorCoord(QMouseEvent*)));
    connect(plotForColocalization->customPlot, SIGNAL(mouseRelease(QMouseEvent*)), this,SLOT(getCursorCoord(QMouseEvent*)));
    connect(plotForActin->customPlot, SIGNAL(mouseRelease(QMouseEvent*)), this,SLOT(getCursorCoord(QMouseEvent*)));

    #ifdef USE_DENS_LATTICES
    plotForDensT = new Plotting(ui->widgetDensT, std::string("Dens free TCR"));;
    plotForDensM = new Plotting(ui->widgetDensM, std::string("Dens free pMHC"));;
    plotForDensTM = new Plotting(ui->widgetDensTM, std::string("Dens TCR-pMHC"));;
    plotForCellKD = new Plotting(ui->widgetCellKD, std::string("CellKD"));;
    connect(plotForDensT->customPlot, SIGNAL(mouseRelease(QMouseEvent*)), this,SLOT(getCursorCoord(QMouseEvent*)));
    connect(plotForDensM->customPlot, SIGNAL(mouseRelease(QMouseEvent*)), this,SLOT(getCursorCoord(QMouseEvent*)));
    connect(plotForDensTM->customPlot, SIGNAL(mouseRelease(QMouseEvent*)), this,SLOT(getCursorCoord(QMouseEvent*)));
    connect(plotForCellKD->customPlot, SIGNAL(mouseRelease(QMouseEvent*)), this,SLOT(getCursorCoord(QMouseEvent*)));
    #endif


    ui->checkBoxICAM1->setChecked(false);
    ui->checkBoxTCR->setChecked(false);
    ui->checkBoxpMHC->setChecked(false);
    ui->checkBoxLFA1->setChecked(false);
    ui->checkBoxTM->setChecked(true);
    ui->checkBoxLI->setChecked(true);
    ui->checkBoxCD2_CD48->setChecked(true);
    ui->checkBoxCD28_CD80->setChecked(true);
    ui->checkBoxPD1_PDL1->setChecked(true);
    ui->checkBoxCD45->setChecked(true);
    ui->checkBoxCD2->setChecked(false);
    ui->checkBoxCD48->setChecked(false);
    ui->checkBoxCD28->setChecked(false);
    ui->checkBoxCD80->setChecked(false);
    ui->checkBoxPD1->setChecked(false);
    ui->checkBoxPDL1->setChecked(false);

    updateColors();
    isRunning = false;
    gamoto = new simulation();
    gamoto->Initialize(true); // is here the problem???
    getDefaultValues();
    t = 0;

    QObject::connect(ui->pushButtonRun, SIGNAL(released()), this, SLOT(runOrStop()));
    QObject::connect(ui->pushButtonReset, SIGNAL(released()), this, SLOT(restart()));
    QObject::connect(ui->pushButtonExportTcell, SIGNAL(released()), this, SLOT(TcellToFig()));
    QObject::connect(ui->pushButtonExportAPC, SIGNAL(released()), this, SLOT(APCToFig()));
    QObject::connect(ui->pushButtonExportColoc, SIGNAL(released()), this, SLOT(ColocalizationToFig()));
    QObject::connect(ui->pushButtonExportTracking, SIGNAL(released()), this, SLOT(TrackingToFig()));
    QObject::connect(ui->pushButtonExportActin, SIGNAL(released()), this, SLOT(ActinToFig()));
    QObject::connect(ui->pushButtonExportALL, SIGNAL(released()), this, SLOT(MakeAllFig()));

    //-------------------------------------------------------------------------------------------------------------------------------
    QObject::connect(ui->checkBoxTCR, SIGNAL(stateChanged(int)), this, SLOT(TCRcheck(int)));
    QObject::connect(ui->checkBoxLFA1, SIGNAL(stateChanged(int)), this, SLOT(LFA1check(int)));
    QObject::connect(ui->checkBoxpMHC, SIGNAL(stateChanged(int)), this, SLOT(pMHCcheck(int)));
    QObject::connect(ui->checkBoxICAM1, SIGNAL(stateChanged(int)), this, SLOT(ICAM1check(int)));
    QObject::connect(ui->checkBoxTM, SIGNAL(stateChanged(int)), this, SLOT(TMcheck(int)));
    QObject::connect(ui->checkBoxCD2_CD48, SIGNAL(stateChanged(int)), this, SLOT(CD2_CD48check(int)));
    QObject::connect(ui->checkBoxCD28_CD80, SIGNAL(stateChanged(int)), this, SLOT(CD28_CD80check(int)));
    QObject::connect(ui->checkBoxPD1_PDL1, SIGNAL(stateChanged(int)), this, SLOT(PD1_PDL1check(int)));
    QObject::connect(ui->checkBoxLI, SIGNAL(stateChanged(int)), this, SLOT(LIcheck(int)));
    QObject::connect(ui->checkBoxCD45, SIGNAL(stateChanged(int)), this, SLOT(CD45check(int)));
    QObject::connect(ui->checkBoxCD2, SIGNAL(stateChanged(int)), this, SLOT(CD2check(int)));
    QObject::connect(ui->checkBoxCD48, SIGNAL(stateChanged(int)), this, SLOT(CD48check(int)));
    QObject::connect(ui->checkBoxCD28, SIGNAL(stateChanged(int)), this, SLOT(CD28check(int)));
    QObject::connect(ui->checkBoxCD80, SIGNAL(stateChanged(int)), this, SLOT(CD80check(int)));
    QObject::connect(ui->checkBoxPD1, SIGNAL(stateChanged(int)), this, SLOT(PD1check(int)));
    QObject::connect(ui->checkBoxPDL1, SIGNAL(stateChanged(int)), this, SLOT(PDL1check(int)));
//-------------------------------------------------------------------------------------------------------------------------------
}


void  MainWindow::getCursorCoord(QMouseEvent *event){
    int x=-1, y = -1;
    if(plotForTcell->parentWidget()->underMouse()){
        x = plotForTcell->customPlot->xAxis->pixelToCoord(event->pos().x());
        y = plotForTcell->customPlot->yAxis->pixelToCoord(event->pos().y());
    }
    if(plotForAPC->parentWidget()->underMouse()){
        x = plotForAPC->customPlot->xAxis->pixelToCoord(event->pos().x());
        y = plotForAPC->customPlot->yAxis->pixelToCoord(event->pos().y());
    }
    if(plotForColocalization->parentWidget()->underMouse()){
        x = plotForColocalization->customPlot->xAxis->pixelToCoord(event->pos().x());
        y = plotForColocalization->customPlot->yAxis->pixelToCoord(event->pos().y());
    }
    if(plotForActin->parentWidget()->underMouse()){
        x = plotForActin->customPlot->xAxis->pixelToCoord(event->pos().x());
        y = plotForActin->customPlot->yAxis->pixelToCoord(event->pos().y());
    }
    #ifdef USE_DENS_LATTICES
    if(plotForDensT->parentWidget()->underMouse()){
        x = plotForDensT->customPlot->xAxis->pixelToCoord(event->pos().x());
        y = plotForDensT->customPlot->yAxis->pixelToCoord(event->pos().y());
    }
    if(plotForDensM->parentWidget()->underMouse()){
        x = plotForDensM->customPlot->xAxis->pixelToCoord(event->pos().x());
        y = plotForDensM->customPlot->yAxis->pixelToCoord(event->pos().y());
    }
    if(plotForDensTM->parentWidget()->underMouse()){
        x = plotForDensTM->customPlot->xAxis->pixelToCoord(event->pos().x());
        y = plotForDensTM->customPlot->yAxis->pixelToCoord(event->pos().y());
    }
    if(plotForCellKD->parentWidget()->underMouse()){
        x = plotForCellKD->customPlot->xAxis->pixelToCoord(event->pos().x());
        y = plotForCellKD->customPlot->yAxis->pixelToCoord(event->pos().y());
    }
    #endif


    //setToolTip(QString("%1 , %2").arg(x).arg(y));
    std::cerr << x << "," << y << std::endl;

    if((x >= 0) && (x < dim) && (y >= 0) && (y < dim)){
        gridpoint* selectedTcell = gamoto->Tcell.access(x,y);
        gridpoint* selectedAPC = gamoto->APC.access(x,y);
        if((!selectedTcell) || (!selectedAPC)) return; //should not happen
        std::cerr << selectedTcell->print() << std::endl;
        std::cerr << selectedAPC->print() << std::endl;
        selectedTcell->isTracked = true;
        selectedAPC->isTracked = true;
    }
}


void MainWindow::getDefaultValues(){

    ui->doubleSpinBoxDifM->setValue(gamoto->Dm);
    ui->doubleSpinBoxDifC->setValue(gamoto->Dc);  // mikro-meters^2 s^(-1)

    ui->doubleSpinBoxWrep->setValue(gamoto->Wrep);
    ui->doubleSpinBoxLrep->setValue(gamoto->Lrep * gamoto->a);
    ui->doubleSpinBoxWatt->setValue(gamoto->Watt);
    ui->doubleSpinBoxLatt->setValue(gamoto->Latt * gamoto->a);

    //ui->doubleSpinBoxKoffTM->setValue(gamoto->Koff_TM );   // s^(-1)    // On- and Off-Probabilities calculation
    //ui->doubleSpinBoxKoffLI->setValue(gamoto->Koff_LI );  // s^(-1)
    //ui->doubleSpinBoxKonTM->setValue(gamoto->Kon_TM ); //2. * pow(10.,4.) * pow(10.,15.); // M^(-1) s^(-1)
    //ui->doubleSpinBoxKonLI->setValue(gamoto->Kon_LI ); // M^(-1) s^(-1)

    ui->spinBoxTCR->setValue(gamoto->init_TCR = 0);
    ui->spinBoxLFA1->setValue(gamoto->init_LFA1 = 0);
    ui->spinBoxpMHC->setValue(gamoto->init_pMHC = 0);
    ui->spinBoxICAM1->setValue(gamoto->init_ICAM = 0);
    ui->spinBoxICAM1->setValue(gamoto->init_CD45 = 0);
    ui->spinBoxICAM1->setValue(gamoto->init_CD2 = 0);
    ui->spinBoxICAM1->setValue(gamoto->init_CD48 = 0);
    ui->spinBoxICAM1->setValue(gamoto->init_CD28 = 0);
    ui->spinBoxICAM1->setValue(gamoto->init_CD80 = 0);

    ui->doubleSpinBoxNucleate->setValue(gamoto->Pnucleate);
    ui->doubleSpinBoxDenucleate->setValue(gamoto->Pdenucleate);
    ui->doubleSpinBoxPolymerize->setValue(gamoto->Ppolymerize);
    ui->doubleSpinBoxDepolymerize->setValue(gamoto->Pdepolymerize);
    ui->doubleSpinBoxRtm->setValue(gamoto->RadiusTM);
    ui->doubleSpinBoxRli->setValue(gamoto->radiusLI);
}

    void MainWindow::restart(){
        /*for (int i=0; i<dim; ++i){
            for (int j=0; j<dim; ++j){
                gamoto->Tcell.CleanLattice(i, j);
                gamoto->APC.CleanLattice(i, j);
            }
        }*/
        t = 0;

        gamoto->Dm = ui->doubleSpinBoxDifM->value();  // mikro-meters^2 s^(-1)
        gamoto->Dc  = ui->doubleSpinBoxDifC->value();  // mikro-meters^2 s^(-1)

        gamoto->Wrep = ui->doubleSpinBoxWrep->value();     // repulsion weight
        gamoto->Lrep = ui->doubleSpinBoxLrep->value()/gamoto->a; //0.4 / a;  // repulsion length /////CHANGEEEEEEEEEEEE SHIT GAMOTO
        gamoto->Watt = ui->doubleSpinBoxWatt->value();    // attraction weight
        gamoto->Latt = ui->doubleSpinBoxLatt->value()/gamoto->a;    // attraction length

        //gamoto->Koff_TM = ui->doubleSpinBoxKoffTM->value();   // s^(-1)    // On- and Off-Probabilities calculation
        //gamoto->Koff_LI = ui->doubleSpinBoxKoffLI->value();  // s^(-1)
        //gamoto->Kon_TM = ui->doubleSpinBoxKonTM->value(); //2. * pow(10.,4.) * pow(10.,15.); // M^(-1) s^(-1)
        //gamoto->Kon_LI = ui->doubleSpinBoxKonLI->value(); // M^(-1) s^(-1)

        gamoto->init_TCR = ui->spinBoxTCR->value();
        gamoto->init_LFA1 = ui->spinBoxLFA1->value();
        gamoto->init_pMHC = ui->spinBoxpMHC->value();
        gamoto->init_ICAM = ui->spinBoxICAM1->value();

        gamoto->Pnucleate = ui->doubleSpinBoxNucleate->value();
        gamoto->Pdenucleate = ui->doubleSpinBoxDenucleate->value();
        gamoto->Ppolymerize = ui->doubleSpinBoxPolymerize->value();
        gamoto->Pdepolymerize = ui->doubleSpinBoxDepolymerize->value();
        gamoto->RadiusTM = ui->doubleSpinBoxRtm->value();
        gamoto->radiusLI = ui->doubleSpinBoxRli->value();

        gamoto->Initialize(false);
    }

    void MainWindow::runOrStop(){
        isRunning = !isRunning;
        ui->pushButtonRun->setText(isRunning ? QString("Stop") : QString("Run"));
        if(isRunning )RunSimulation();
    }

    void MainWindow::TcellToFig(QString fileNameTcell){
        if(fileNameTcell.size() < 1) fileNameTcell = QFileDialog::getSaveFileName(this, QString("File to save image"), QString("~"), QString("*.png"));
        plotForTcell->customPlot->savePng(fileNameTcell);
    }
    void MainWindow::APCToFig(QString fileNameAPC){
        if(fileNameAPC.size() < 1) fileNameAPC = QFileDialog::getSaveFileName(this, QString("File to save image"), QString("~"), QString("*.png"));
        plotForAPC->customPlot->savePng(fileNameAPC);
    }
    void MainWindow::ColocalizationToFig(QString fileNameColocalization){
        if(fileNameColocalization.size() < 1) fileNameColocalization = QFileDialog::getSaveFileName(this, QString("File to save image"), QString("~"), QString("*.png"));
        plotForColocalization->customPlot->savePng(fileNameColocalization);
    }
    void MainWindow::TrackingToFig(QString fileNameTracking){
        if(fileNameTracking.size() < 1) fileNameTracking = QFileDialog::getSaveFileName(this, QString("File to save image"), QString("~"), QString("*.png"));
        plotForTracking->customPlot->savePng(fileNameTracking);
    }
    void MainWindow::ActinToFig(QString fileNameActin){
        if(fileNameActin.size() < 1) fileNameActin = QFileDialog::getSaveFileName(this, QString("File to save image"), QString("~"), QString("*.png"));
        plotForActin->customPlot->savePng(fileNameActin);
    }

    #ifdef USE_DENS_LATTICES
    void MainWindow::DensTToFig(QString fileNameDensT){
        if(fileNameDensT.size() < 1) fileNameDensT = QFileDialog::getSaveFileName(this, QString("File to save image"), QString("~"), QString("*.png"));
        plotForDensT->customPlot->savePng(fileNameDensT);
    }
    void MainWindow::DensMToFig(QString fileNameDensM){
        if(fileNameDensM.size() < 1) fileNameDensM = QFileDialog::getSaveFileName(this, QString("File to save image"), QString("~"), QString("*.png"));
        plotForDensM->customPlot->savePng(fileNameDensM);
    }
    void MainWindow::DensTMToFig(QString fileNameDensTM){
        if(fileNameDensTM.size() < 1) fileNameDensTM = QFileDialog::getSaveFileName(this, QString("File to save image"), QString("~"), QString("*.png"));
        plotForDensTM->customPlot->savePng(fileNameDensTM);
    }
    void MainWindow::CellKDToFig(QString fileNameCellKD){
        if(fileNameCellKD.size() < 1) fileNameCellKD = QFileDialog::getSaveFileName(this, QString("File to save image"), QString("~"), QString("*.png"));
        plotForCellKD->customPlot->savePng(fileNameCellKD);
    }
    #endif

        void MainWindow::MakeAllFig(){
            QString suffix = QString::number(gamoto->T);
            QString folderToSave = QFileDialog::getExistingDirectory( this, QString("File to save image"), QString("~"));
            TcellToFig(folderToSave + QString("Tcell") + suffix + QString(".png"));
            APCToFig(folderToSave + QString("APC") + suffix + QString(".png"));
            ColocalizationToFig(folderToSave + QString("Coloc") + suffix + QString(".png"));
            TrackingToFig(folderToSave + QString("Tracking") + suffix + QString(".png"));
            ActinToFig(folderToSave + QString("Actin") + suffix + QString(".png"));
            #ifdef USE_DENS_LATTICES
            DensTToFig(folderToSave + QString("DensT") + suffix + QString(".png"));
            DensMToFig(folderToSave + QString("DensM") + suffix + QString(".png"));
            DensTMToFig(folderToSave + QString("DensTM") + suffix + QString(".png"));
            CellKDToFig(folderToSave + QString("CellKD") + suffix + QString(".png"));
            #endif
        }


//-------------------------------------------------------------------------------------------------------------------------------
    void MainWindow::updateColors(){
        QMap<double, QColor> newColor = AdaptColorScale(ui->checkBoxTCR->isChecked(), ui->checkBoxLFA1->isChecked(), ui->checkBoxpMHC->isChecked(), ui->checkBoxICAM1->isChecked(),
                                                        ui->checkBoxTM->isChecked(), ui->checkBoxLI->isChecked(), ui->checkBoxCD45->isChecked(),
                                                        ui->checkBoxCD2->isChecked(), ui->checkBoxCD48->isChecked(), ui->checkBoxCD2_CD48->isChecked(),
                                                        ui->checkBoxCD28->isChecked(), ui->checkBoxCD80->isChecked(), ui->checkBoxCD28_CD80->isChecked(),
                                                        ui->checkBoxPD1->isChecked(), ui->checkBoxPDL1->isChecked(), ui->checkBoxPD1_PDL1->isChecked());
        plotForTcell->ChangeColors(newColor);
        plotForAPC->ChangeColors(newColor);
        plotForColocalization->ChangeColors(newColor);
        plotForTracking->ChangeColors(newColor);
        plotForActin->ChangeColors(newColor);

        #ifdef USE_DENS_LATTICES
        QMap<double, QColor> colStops;
        //#define coefHere 10./27.0
        colStops.insert(0.0,  QColor(Qt::black)); // nothing
        colStops.insert(0.002,  QColor(Qt::blue)); // nothing
        colStops.insert(0.05,  QColor(Qt::cyan)); // nothing
        colStops.insert(0.1,  QColor(Qt::yellow));// border - means 1
        colStops.insert(1.,  QColor(Qt::red)); // means 10
        colStops.insert(100.,  QColor(Qt::green)); // out of bounds

//        colStops.insert(0.0,  QColor(Qt::black)); // nothing
//        colStops.insert(0.0025,  QColor(Qt::darkBlue)); // nothing
//        colStops.insert(0.005,  QColor(Qt::blue)); // nothing
//        colStops.insert(0.0075,  QColor(Qt::darkCyan)); // nothing
//        colStops.insert(0.01,  QColor(Qt::cyan)); // nothing
//        colStops.insert(0.025,  QColor(Qt::darkYellow)); // nothing
//        colStops.insert(0.05,  QColor(Qt::yellow)); // nothing
//        colStops.insert(0.075,  QColor(Qt::red)); // nothing
//        colStops.insert(0.1,  QColor(Qt::darkRed));// border - means 1
//        colStops.insert(1.,  QColor(Qt::darkRed)); // means 10
//        colStops.insert(1000.,  QColor(Qt::green)); // out of bounds
        plotForDensT->ChangeColors(colStops);
        plotForDensT->colorScale->setDataRange(QCPRange(0.0, 10));
        plotForDensM->ChangeColors(colStops);
        plotForDensM->colorScale->setDataRange(QCPRange(0.0, 10));
        plotForDensTM->ChangeColors(colStops);
        plotForDensTM->colorScale->setDataRange(QCPRange(0.0, 10));
        plotForCellKD->ChangeColors(colStops);
        plotForCellKD->colorScale->setDataRange(QCPRange(0.0, 10));
        #endif
    }

    void MainWindow::TCRcheck(int ){
        updateColors();
    }
    void MainWindow::LFA1check(int){
        updateColors();
    }
    void MainWindow::pMHCcheck(int ){
        updateColors();
    }
    void MainWindow::ICAM1check(int){
        updateColors();
    }
    void MainWindow::TMcheck(int){
        updateColors();
    }
    void MainWindow::LIcheck(int){
        updateColors();
    }
    void MainWindow::CD2_CD48check(int){
        updateColors();
    }
    void MainWindow::CD28_CD80check(int){
        updateColors();
    }
    void MainWindow::PD1_PDL1check(int){
        updateColors();
    }
    void MainWindow::CD45check(int){
        updateColors();
    }
    void MainWindow::CD2check(int){
        updateColors();
    }
    void MainWindow::CD48check(int){
        updateColors(); 
    }
    void MainWindow::CD28check(int){
        updateColors();
    }
    void MainWindow::CD80check(int){
        updateColors();
    }
    void MainWindow::PD1check(int){
        updateColors();
    }
    void MainWindow::PDL1check(int){
        updateColors();
    }
//-------------------------------------------------------------------------------------------------------------------------------
    void MainWindow::Koff_LI(double val1){
        ui->doubleSpinBoxKoffLI->setValue(val1);
    }
    void MainWindow::Koff_TM(double val2){
        ui->doubleSpinBoxKoffTM->setValue(val2);
    }
    void MainWindow::Kon_LI(double val3){
        ui->doubleSpinBoxKonLI->setValue(val3);
    }
    void MainWindow::Kon_TM(double val4){
        ui->doubleSpinBoxKonTM->setValue(val4);
    }
//-------------------------------------------------------------------------------------------------------------------------------
    void MainWindow::RunSimulation(){
        //new simulation starting

        restart();

        QString fileFolder = QFileDialog::getExistingDirectory(this, QString("Folder to save pictures"), QString("~"));
        std::ofstream vid(fileFolder.toStdString() + "/videoT_APC.sh");
        vid << "ffmpeg -framerate 100  -start_number 10001 -i \"Tcell%05d.png\" -c:v libx264 -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" -r 30  -pix_fmt yuv420p Tcellout.avi\n";
        vid << "ffmpeg -framerate 100  -start_number 10001 -i \"APC%05d.png\" -c:v libx264 -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" -r 30  -pix_fmt yuv420p APCout.avi\n";
        vid << "ffmpeg -framerate 100  -start_number 10001 -i \"Actin%05d.png\" -c:v libx264 -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" -r 30  -pix_fmt yuv420p Actinout.avi\n";
        vid << "ffmpeg -framerate 100  -start_number 10001 -i \"Tracking%05d.png\" -c:v libx264 -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" -r 30  -pix_fmt yuv420p Trackingout.avi\n";
        vid << "ffmpeg -framerate 100  -start_number 10001 -i \"Colocalization%05d.png\" -c:v libx264 -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" -r 30  -pix_fmt yuv420p Colocalizationout.avi\n";

        vid.close();
        std::string commandOpt = std::string("chmod +x \"") + fileFolder.toStdString() + "/videoT_APC.sh\"";
        std::cerr << "trying " << commandOpt << std::endl;
        system(commandOpt.c_str());

        std::cout << "saving in folder " << fileFolder.toStdString() << std::endl;

        TimerStart();
        int count = 0;
        double Showtime = 1.;
        double StopTime = 600.;
        while(isRunning){
            //std::cerr << TimeEvolution * 0.5 << endl;
            gamoto-> Simulate(Showtime);
            count++;
            /*if(!(count % 100)){
                std::cout<< "Simulation time "<< t << " Time passed = "<< TimerStop()<<std::endl;
                //exit(-1);
            }*/
            if(t > StopTime){
                std::cout<< "Time passed = "<< TimerStop()<<std::endl;
                exit(-1);
            }
            t += Showtime;
            TcellToFig(fileFolder + QString("/Tcell") + QString::number((int)(10000+t)) + QString(".png"));
            APCToFig(fileFolder + QString("/APC") +  QString::number((int)(10000+t)) + QString(".png"));
            ColocalizationToFig(fileFolder + QString("/Colocalization") + QString::number(10000+t) + QString(".png"));
            TrackingToFig(fileFolder + QString("/Tracking") + QString::number((int)(10000+t)) + QString(".png"));

            plotForTcell->PlotLattice(&gamoto->Tcell);
            plotForAPC->PlotLattice(&gamoto->APC);
            plotForColocalization->PlotLattice(&gamoto->Colocalization);
//            plotForActin->PlotLattice(&gamoto->Actin);
            plotForTracking->PlotLattice(&gamoto->Tracker);

//            #ifdef USE_DENS_LATTICES
//            plotForDensT->PlotLattice(&gamoto->DensT);
//            plotForDensM->PlotLattice(&gamoto->DensM);
//            plotForDensTM->PlotLattice(&gamoto->DensTM);
//            plotForCellKD->PlotLattice(&gamoto->CellKD);
//            #endif

            ui->lcdNumber->display(t);
            QCoreApplication::processEvents();
        }
    }

MainWindow::~MainWindow()
{
    delete ui;
}

    Plotting::Plotting(QWidget *parent, std::string title) :
        QWidget(parent){
        customPlot = new QCustomPlot(parent);
        customPlot->resize(parent->size());

//        customPlot->plotLayout()->insertRow(0);
       // customPlot->plotLayout()->insertColumn(0);
        customPlot->plotLayout()->addElement(1,0, new QCPPlotTitle(customPlot,title.c_str()));

        // configure axis rect:
        customPlot->setInteractions(QCP::iRangeDrag|QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
        customPlot->axisRect()->setupFullAxesBox(true);
        //customPlot->xAxis->setLabel("x");
        //customPlot->yAxis->setLabel("y");

        // set up the QCPColorMap:
        colorMap = new QCPColorMap(customPlot->xAxis, customPlot->yAxis);
        //(phi)colorMap->setDataRange(QCPRange(0.0, 8.0));

        customPlot->addPlottable(colorMap);
        int nx = dim;
        int ny = dim;

        colorMap->setInterpolate(true);///////// To decide smooth or not pixels
        colorMap->data()->setSize(nx, ny); // we want the color map to have nx * ny data points
        colorMap->data()->setRange(QCPRange(0, dim), QCPRange(0, dim)); // and span the coordinate range in both key (x) and value (y) dimensions

        colorScale = new QCPColorScale(customPlot);
        colorScale->setDataRange(QCPRange(0.0, 27.0));      // apparently, the gradient is always between 0 and 1, and this does the scale
        colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
        colorScale->axis()->setLabel("Molecule Type");

        customPlot->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
        colorMap->setColorScale(colorScale); // associate the color map with the color scale

        //myGrad->setLevelCount(21);                  // number of color levels in the gradient
        QMap<double, QColor> colStops;
        colStops.insert(0.0, QColor(Qt::black));      // nothing
        colStops.insert(0.05, QColor(Qt::black));     // border
        colStops.insert(0.1, QColor(Qt::black));
        colStops.insert(0.15, QColor(Qt::black));     // TCR
        colStops.insert(0.2, QColor(Qt::black));
        colStops.insert(0.25, QColor(Qt::black));     // LFA1
        colStops.insert(0.3, QColor(Qt::black));
        colStops.insert(0.35, QColor(Qt::black));     // pMHC
        colStops.insert(0.4, QColor(Qt::black));
        colStops.insert(0.45, QColor(Qt::black));     // ICAM1
        colStops.insert(0.5, QColor(Qt::black));
        colStops.insert(0.55, QColor(Qt::darkGreen)); // TCR-pMHC
        colStops.insert(0.6, QColor(Qt::darkGreen));
        colStops.insert(0.65, QColor(Qt::red));       // LFA1-ICAM1
        colStops.insert(0.7, QColor(Qt::red));
        colStops.insert(0.75, QColor(Qt::yellow));    // Nucleation Point
        colStops.insert(0.8, QColor(Qt::yellow));
        colStops.insert(0.85, QColor(Qt::magenta));   // NucleationPolymerized Point
        colStops.insert(0.9, QColor(Qt::magenta));
        colStops.insert(0.95, QColor(Qt::green));     // Polymerization Point
        colStops.insert(1.0, QColor(Qt::green));
        colStops.insert(1.05, QColor(Qt::white));     //cd45
        colStops.insert(1.1, QColor(Qt::white));
        colStops.insert(1.15, QColor(Qt::black));     //empty actin
        colStops.insert(1.2, QColor(Qt::black));
        colStops.insert(1.25, QColor(Qt::cyan));      //cd2
        colStops.insert(1.3, QColor(Qt::cyan));
        colStops.insert(1.35, QColor(Qt::yellow));    //cd48
        colStops.insert(1.4, QColor(Qt::yellow));
        colStops.insert(1.45, QColor(Qt::yellow));    //cd2-cd48
        colStops.insert(1.5, QColor(Qt::yellow));

        ChangeColors(colStops);

    }

    void Plotting::PlotLattice(doubleLattice *Cell){
        int nx = dim;
        int ny = dim;
        //colorMap->data()->clear();

        // now we assign some data, by accessing the QCPColorMapData instance of the color map:
        for (int xIndex=0; xIndex<nx; ++xIndex){
          for (int yIndex=0; yIndex<ny; ++yIndex){
            colorMap->data()->setCell(xIndex, yIndex, Cell->access(xIndex,yIndex));
          }
        }
        customPlot->rescaleAxes();   // rescale the key (x) and value (y) axes so the whole color map is visible:
        customPlot->replot();
    }

    void Plotting::PlotLattice(lattice *Cell){
        int nx = dim;
        int ny = dim;
        //colorMap->data()->clear();

        // now we assign some data, by accessing the QCPColorMapData instance of the color map:
        double z;
        for (int xIndex=0; xIndex<nx; ++xIndex){
          for (int yIndex=0; yIndex<ny; ++yIndex){
            //colorMap->data()->cellToCoord(xIndex, yIndex, &x, &y);
            //double r = 3*qSqrt(x*x+y*y)+1e-2;
            z = Cell->access(xIndex,yIndex)->getAgent(); // (phi)
            //std::cout << z << std::endl;
            colorMap->data()->setCell(xIndex, yIndex, z);
          }
        }

        // rescale the key (x) and value (y) axes so the whole color map is visible:
        customPlot->rescaleAxes();

        customPlot->replot();
    }

    void Plotting::ChangeColors(QMap<double, QColor> colStops){

        myGrad = new QCPColorGradient();      // when doesn't say a gradient, he chooses gpGold

//-------------------------------------------------------------------------------------------------------------------------------------
        myGrad->setLevelCount(51);       //CHANGE HERE WHEN YOU ADD MOLECULES           // number of color levels in the gradient
//-------------------------------------------------------------------------------------------------------------------------------------
        myGrad->setColorStops(colStops);

        colorScale->setGradient(*myGrad);
        // set the color gradient of the color map to one of the presets:
        // colorMap->setGradient(QCPColorGradient::gpPolar);
        // we could have also created a QCPColorGradient instance and added own colors to
        // the gradient, see the documentation of QCPColorGradient for what's possible.

    }


#ifdef UseGraphInterface

int main(int argc, char** argv){
    QApplication skata(argc, argv);
    MainWindow* mw = new MainWindow();
    mw->show();

    //lattice* Tcell = new lattice(); // create the lattice for the T cell
    //lattice* APC = new lattice();   // create the lattice for the APC

    //lattice Tcell;
    //lattice APC;
    //main_2(*Tcell, *APC);

    return skata.exec();
}
#endif


