#ifndef POSTPROCESSING_HH
#define POSTPROCESSING_HH

#include <string>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <filesystem>

#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TLeafC.h>
#include <TH1.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TError.h>

#include "Configuration.hh"
#include "GenSurface.hh"

class TFile;
class TH1;

class PostProcessing {
public:
    PostProcessing(std::string outputFolderName,
                   std::string particle,
                   double norm_cm2,
                   bool isotropic);

    ~PostProcessing();

    void ExtractNtData();
    void SaveResponse();

private:
    std::string outputFolderName;
    double eMinMeV;
    double eMaxMeV;
    std::string particleName;
    double norm;
    bool isotropic;

    std::unique_ptr<TFile> rootFile;

    std::string postProcessingDir;
    std::string runDir;
    std::string csvDir;
    std::string histogramsDir;

    void OpenRootFile();
    void PrepareOutputDirs();

    void ExportTreeToCsv(const std::string& treeName,
                         const std::string& csvPath);

    TH1* GetHist(const std::string& histName);
    TH1* GetHistOrThrow(const std::string& histName);

    void SaveHistPng(const std::string& histName,
                     const std::string& outPngPath,
                     const std::string& plotTitle,
                     const std::string& yTitle,
                     bool filled,
                     bool log_y = false,
                     bool weighted = false);

    [[nodiscard]] std::string AreaUnit() const;
    [[nodiscard]] std::string AreaUnitRoot() const;

    static double BinCentre(double eLow, double eHigh);
    static double EffAreaErrFromCounts(double n0, double n, double effArea);
};


#endif //POSTPROCESSING_HH
