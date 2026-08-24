#include "Loader.hh"

using namespace Configuration;

Loader::Loader(int argc, char** argv) {
    numThreads = G4Threading::G4GetNumberOfCores();
    useUI = true;
    macroFile = "../run.mac";
    detectorType = "CsI";
    fluxType = "Uniform";
    fluxDirection = "isotropic";
    eThreshold = 0 * MeV;
    outputFile = "GeigerCounter.root";
    nBins = 1000;
    saveSecondaries = false;

    auto value = [argc, argv](const int i) -> std::string {
        if (i + 1 >= argc) {
            G4Exception("Loader::Loader", "MISSING_ARGUMENT", FatalException,
                        (std::string("Option ") + argv[i] + " requires a value").c_str());
            return "";
        }
        return argv[i + 1];
    };

    for (int i = 1; i < argc; i++) {
        const std::string input(argv[i]);

        if (input == "-i" || input == "--input") {
            macroFile = value(i);
            useUI = false;
        } else if (input == "-t" || input == "--threads") {
            numThreads = std::stoi(value(i));
        } else if (input == "--bins") {
            nBins = std::stoi(value(i));
        } else if (input == "-noUI") {
            useUI = false;
        } else if (input == "-d" || input == "--detector") {
            detectorType = value(i);
        } else if (input == "-f" || input == "--flux-type") {
            fluxType = value(i);
        } else if (input == "--flux-dir" || input == "--f-dir" || input == "-fd") {
            fluxDirection = value(i);
        } else if (input == "--save-secondaries") {
            saveSecondaries = true;
        } else if (input == "-o" || input == "--output-file") {
            outputFile = value(i);
            outputFile += ".root";
        }
    }

    configPath = "../Flux_config/" + fluxType + "_params.txt";

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(time(nullptr));

#ifdef G4MULTITHREADED
    runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(numThreads);
#else
    runManager = new G4RunManager;
#endif

    auto* realWorld = new Geometry();
    runManager->SetUserInitialization(realWorld);
    auto* physicsList = new QBBC;
    physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
    physicsList->ReplacePhysics(new G4RadioactiveDecayPhysics());

    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);

    Emin = std::max({std::stod(ReadValue("E_min:", "")) * MeV, eThreshold});
    Emax = std::stod(ReadValue("E_max:", "")) * MeV;
    if (fluxType == "Table") {
        std::string path = ReadValue("table_path:", "");
        auto energyTable = Utils::ReadCSV(path, 1., false, MeV);

        Emin = std::max({energyTable.GetMinE(), Emin});
        Emax = std::min({energyTable.GetMaxE(), Emax});
    }

    if (fluxType == "Uniform") {
        std::string isLogStr = ReadValue("is_log:", configPath);
        isLogBin = isLogStr == "1" || isLogStr == "true";
    }

    genSurface = GenSurface::For(fluxDirection);
    area = genSurface.Norm_cm2();
    runManager->SetUserInitialization(new ActionInitialization(area));
    runManager->Initialize();

    visManager = new G4VisExecutive;
    visManager->Initialize();
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (!useUI) {
        const G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macroFile);
    } else {
        auto* ui = new G4UIExecutive(argc, argv, "qt");
        UImanager->ApplyCommand("/control/execute ../vis.mac");
        ui->SessionStart();
        delete ui;
    }

    const auto* runAction = dynamic_cast<const RunAction*>(runManager->GetUserRunAction());
    if (runAction) {
        const auto& [t1Only, t2Only, t1AndT2] = runAction->GetCounts();
        tube1Only = t1Only;
        tube2Only = t2Only;
        both = t1AndT2;
        nGenerated = runAction->GetNGenerated();
        effArea1 = runAction->GetEffArea();
        effAreaTel = runAction->GetEffAreaTelescope();
    }
    SaveConfig();
    RunPostProcessing();
}

Loader::~Loader() {
    delete runManager;
    delete visManager;
}


std::string Loader::ReadValue(const std::string& key, const std::string& filepath = "") const {
    std::ifstream file(filepath.empty() ? configPath : filepath);
    if (!file.is_open()) {
        G4Exception("Loader::ReadValue", "FILE_OPEN_FAIL",
                    FatalException, ("Cannot open " + (filepath.empty() ? configPath : filepath)).c_str());
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.find(key) != std::string::npos) {
            return line.substr(key.length() + 1);
        }
    }

    return "";
}


inline std::string Trim(std::string st) {
    auto notSpace = [](const unsigned char c) {
        return !std::isspace(c);
    };
    st.erase(st.begin(), std::find_if(st.begin(), st.end(), notSpace));
    st.erase(std::find_if(st.rbegin(), st.rend(), notSpace).base(), st.end());
    return st;
}


std::vector<G4String> Split(const G4String& line) {
    std::vector<G4String> result;
    std::stringstream ss(line);
    G4String token;
    while (std::getline(ss, token, ',')) {
        token = Trim(token);
        if (!token.empty())
            result.push_back(token);
    }
    return result;
}


void Loader::SaveConfig() const {
    const int N = nGenerated;

    EnergyRange er{};
    FluxType fType{};
    FluxParams fp{};

    er.Emin = Emin;
    er.Emax = Emax;

    if (fluxType == "PLAW") {
        fType = FluxType::PLAW;
        fp.A = std::stod(ReadValue("A:"));
        fp.alpha = std::stod(ReadValue("alpha:"));
        fp.E_piv = std::stod(ReadValue("E_Piv:"));
    } else if (fluxType == "COMP") {
        fType = FluxType::COMP;
        fp.A = std::stod(ReadValue("A:"));
        fp.alpha = std::stod(ReadValue("alpha:"));
        fp.E_piv = std::stod(ReadValue("E_Piv:"));
        fp.E_peak = std::stod(ReadValue("E_Peak:"));
    } else if (fluxType == "SEP") {
        fType = FluxType::SEP;
        fp.sep_year = std::stoi(ReadValue("year:"));
        fp.sep_order = std::stoi(ReadValue("order:"));
        fp.sep_csv_path = "../SEP_coefficients.CSV";
    } else if (fluxType == "Galactic") {
        fType = FluxType::GALACTIC;
        fp.phiMV = std::stod(ReadValue("phiMV:"));
        fp.particle = ReadValue("particle:");
    } else if (fluxType == "Table") {
        fType = FluxType::TABLE;
        fp.particle = ReadValue("particle:");
        fp.table_path = ReadValue("table_path:");
    } else {
        fType = FluxType::UNIFORM;
    }

    RateCounts counts{tube1Only, tube2Only, both};

    RateResult rr{};
    bool rate_ok = true;
    try {
        rr = computeRate(fType, fp, er, area, N, counts);
    }
    catch (const std::exception& ex) {
        rate_ok = false;
    }

    RateResult rrReal{};
    bool rate_real_ok = true;
    try {
        rrReal = computeRateReal(fType, fp, er, effArea1, effAreaTel, nBins);
    }
    catch (const std::exception& ex) {
        rate_real_ok = false;
    }

    std::ostringstream buf;

    buf << "N_generated: " << N << "\n";
    buf << "N_bins: " << nBins << "\n";
    buf << "Log_binning: " << (isLogBin ? 1 : 0) << "\n\n";
    buf << "Flux_type: " << fluxType << "\n";
    buf << "Flux_dir: " << fluxDirection << "\n";

    buf << "Flux_params:\n{\n\t";
    if (fluxType == "PLAW") {
        buf << "A: " << std::stod(ReadValue("A:")) << ",\n\t";
        buf << "alpha: " << std::stod(ReadValue("alpha:")) << ",\n\t";
        buf << "E_Piv: " << std::stod(ReadValue("E_Piv:")) << " MeV\n";
    } else if (fluxType == "COMP") {
        buf << "A: " << std::stod(ReadValue("A:")) << ",\n\t";
        buf << "alpha: " << std::stod(ReadValue("alpha:")) << ",\n\t";
        buf << "E_Piv: " << std::stod(ReadValue("E_Piv:")) << " MeV,\n\t";
        buf << "E_Peak: " << std::stod(ReadValue("E_Peak:")) << " MeV\n";
    } else if (fluxType == "SEP") {
        buf << "year: " << std::stoi(ReadValue("year:")) << ",\n\t";
        buf << "order: " << std::stoi(ReadValue("order:")) << "\n";
    } else if (fluxType == "Galactic") {
        buf << "phiMV: " << std::stod(ReadValue("phiMV:")) << " MV,\n\t";
        buf << "particle: " << ReadValue("particle:") << "\n";
    } else if (fluxType == "Table") {
        buf << "table_path: " << ReadValue("table_path:") << ",\n\t";
        buf << "particle: " << ReadValue("particle:") << "\n";
    } else if (fluxType == "Uniform") {
        buf << "fractions: " << ReadValue("fractions:") << "\n";
    }
    buf << "}\n\n";

    buf << "Particles: [";
    if (fluxType == "PLAW") {
        buf << "gamma";
    } else if (fluxType == "SEP") {
        buf << "proton";
    } else if (fluxType == "Galactic" or fluxType == "Table") {
        buf << ReadValue("particle:");
    } else if (fluxType == "Uniform") {
        buf << ReadValue("particles:");
    }
    buf << "]\n";

    buf << "Energies:\n{\n\t";
    if (fluxType == "PLAW" || fluxType == "COMP") {
        buf << "gamma: ";
    } else if (fluxType == "SEP") {
        buf << "proton: ";
    } else if (fluxType == "Galactic" or fluxType == "Table") {
        buf << ReadValue("particle:") << ": ";
    } else if (fluxType == "Uniform") {
        std::vector<G4String> particles = Split(ReadValue("particles:"));
        std::vector<G4String> EminVec = Split(ReadValue("E_min:"));
        std::vector<G4String> EmaxVec = Split(ReadValue("E_max:"));
        for (size_t i = 0; i < particles.size(); i++) {
            buf << (i == 0 ? "" : "\t") << particles[i] << ": (" << EminVec[i] << " MeV, " << EmaxVec[i] << " MeV),\n";
        }
    }
    if (fluxType != "Uniform") {
        buf << "(" << Emin << " MeV, " << Emax << " MeV)\n";
    }
    buf << "}\n\n";

    const bool isotropic = genSurface.IsIsotropic();
    const std::string area_dim = isotropic ? " cm^2*sr" : " cm^2";
    const std::string area_dim_inv = isotropic ? " cm^-2*sr^-1" : " cm^-2";

    buf << "Detector:\n{\n\t";
    buf << "gas_pressure: " << gasPressure / (atmosphere / 760.) << " torr,\n\t";
    buf << "gas_temperature: " << gasTemperature / kelvin << " K,\n\t";
    buf << "range_cut: " << counterRangeCut / um << " um,\n\t";
    buf << "threshold: " << eThreshold / eV << " eV\n}\n\n";

    buf << "Generation_surface:\n{\n\t";
    buf << "shape: " << genSurface.ShapeName() << ",\n\t";
    if (isotropic) {
        buf << "R_gen: " << genSurface.Radius() / mm << " mm,\n\t";
    } else {
        buf << "half_u: " << genSurface.HalfU() / mm << " mm,\n\t";
        buf << "half_v: " << genSurface.HalfV() / mm << " mm,\n\t";
        buf << "standoff: " << genSurface.Standoff() / mm << " mm,\n\t";
        buf << "axis: (" << genSurface.Axis().x() << ", " << genSurface.Axis().y()
            << ", " << genSurface.Axis().z() << "),\n\t";
    }
    buf << "S_perp: " << genSurface.SPerp_cm2() << " cm^2,\n\t";
    buf << "geometric_factor: " << genSurface.GeomFactor_cm2sr() << " cm^2*sr,\n\t";
    buf << "normalisation: " << area << area_dim << "\n}\n\n";

    buf << "Counts:\n{\n\t";
    buf << "N_1: " << tube1Only + both << ",\n\t";
    buf << "N_tel: " << both << ",\n\t";
    buf << "tube1_only: " << tube1Only << ",\n\t";
    buf << "tube2_only: " << tube2Only << ",\n\t";
    buf << "both: " << both << "\n}\n\n";

    buf << std::fixed << std::setprecision(6);

    buf << "Response:\n{\n\t";
    if (N > 0) {
        buf << "A_1: " << area * (tube1Only + both) / static_cast<double>(N) << area_dim << ",\n\t";
        buf << "A_tel: " << area * both / static_cast<double>(N) << area_dim << "\n}\n\n";
    } else {
        buf << "A_1: NaN,\n\tA_tel: NaN\n}\n\n";
    }

    buf << "Rates:\n{\n\t";
    if (rate_ok) {
        buf << "Integral: " << rr.integral << area_dim_inv << " * s^-1,\n\t";
        buf << "Ndot: " << rr.Ndot << " s^-1,\n\t";
        buf << "Rate_1: " << rr.rate1 << " s^-1,\n\t";
        buf << "Rate_tel: " << rr.rateTel << " s^-1,\n\t";
    } else {
        buf << "Integral: NaN,\n\tNdot: NaN,\n\tRate_1: NaN,\n\tRate_tel: NaN,\n\t";
    }
    if (rate_real_ok) {
        buf << "Rate_real_1: " << rrReal.rateReal1 << " s^-1,\n\t";
        buf << "Rate_real_tel: " << rrReal.rateRealTel << " s^-1\n";
    } else {
        buf << "Rate_real_1: NaN,\n\tRate_real_tel: NaN\n";
    }
    buf << "}\n\n";

    auto sanitize = [](std::string ss) {
        for (char& c : ss) if (c == ' ') c = '_';
        return ss;
    };

    std::string filename = "info_" + detectorType + "_" + fluxType;
    if (fluxType == "Galactic") {
        const std::string part = ReadValue("particle:");
        const std::string phi = ReadValue("phiMV:");
        filename += "_particle:" + part + "_phiMV:" + phi + ".txt";
    } else if (fluxType == "Uniform") {
        const std::string part = ReadValue("particles:");
        filename += "_particle:" + part + ".txt";
    } else {
        filename += ".txt";
    }
    filename = sanitize(filename);

    std::ofstream out(filename);
    if (!out.is_open()) {
        G4cerr << "Ошибка: не удалось открыть файл " << filename << G4endl;
        return;
    }
    out << buf.str();
    out.close();

    std::cout << "Configuration saved in " << filename << std::endl;
}


void Loader::RunPostProcessing() const {
    auto sanitize = [](std::string ss) {
        for (char& c : ss) if (c == ' ') c = '_';
        return ss;
    };
    std::string part;

    try {
        std::cout << "Processing... ";
        std::string outDir = fluxType;
        if (fluxType == "Galactic") {
            const std::string phi = ReadValue("phiMV:");
            part = ReadValue("particle:");
            outDir += "_" + part + "_phiMV:" + phi;
        } else if (fluxType == "Uniform") {
            part = ReadValue("particles:");
            outDir += "_" + part;
        } else if (fluxType == "PLAW" || fluxType == "COMP") {
            part = "gamma";
            outDir += "_" + part;
        } else if (fluxType == "SEP") {
            part = "proton";
            outDir += "_" + part;
        } else if (fluxType == "Table") {
            part = ReadValue("particle:");
            outDir += "_" + part;
        }
        outDir += "_" + fluxDirection;
        outDir = sanitize(outDir);
        part = sanitize(part);

        PostProcessing postProcessing(outDir, part, area, genSurface.IsIsotropic());

        postProcessing.ExtractNtData();
        if (Emin < Emax) {
            postProcessing.SaveResponse();
        }

        std::cout << "Done!\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}
