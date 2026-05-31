#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct Accumulator {
    string vtMean;
    string vtStd;
    string vrMean;
    string vrStd;
    double sumA2 = 0.0;
    double sumM2LogScore = 0.0;
    int count = 0;
};

static vector<string> splitCsvLine(const string& line) {
    vector<string> fields;
    string current;
    bool inQuotes = false;

    for (size_t i = 0; i < line.size(); i++) {
        char c = line[i];
        if (c == '"') {
            if (inQuotes && i + 1 < line.size() && line[i + 1] == '"') {
                current += '"';
                i++;
            } else {
                inQuotes = !inQuotes;
            }
        } else if (c == ',' && !inQuotes) {
            fields.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }

    fields.push_back(current);
    return fields;
}

static int columnIndex(const vector<string>& header, const string& name) {
    for (int i = 0; i < static_cast<int>(header.size()); i++) {
        if (header[i] == name) {
            return i;
        }
    }

    return -1;
}

static bool hasColumnIndexes(const vector<int>& indexes) {
    for (int index : indexes) {
        if (index < 0) {
            return false;
        }
    }

    return true;
}

static string keyFor(
    const string& vtMean,
    const string& vtStd,
    const string& vrMean,
    const string& vrStd
) {
    return vtMean + '\t' + vtStd + '\t' + vrMean + '\t' + vrStd;
}

static void mergeOneFile(const string& inputPath, const string& outputPath) {
    ifstream fin(inputPath);
    if (!fin) {
        cerr << "Missing input: " << inputPath << endl;
        ofstream emptyOutput(outputPath);
        emptyOutput << "vt_mean,vt_std,vr_mean,vr_std,mean_A2,mean_m2_log_score,merged_count\n";
        return;
    }

    string line;
    if (!getline(fin, line)) {
        cerr << "Empty input: " << inputPath << endl;
        return;
    }

    vector<string> header = splitCsvLine(line);
    int vtMeanIndex = columnIndex(header, "vt_mean");
    int vtStdIndex = columnIndex(header, "vt_std");
    int vrMeanIndex = columnIndex(header, "vr_mean");
    int vrStdIndex = columnIndex(header, "vr_std");
    int meanA2Index = columnIndex(header, "mean_A2");
    int meanM2Index = columnIndex(header, "mean_m2_log_score");

    if (!hasColumnIndexes({
            vtMeanIndex,
            vtStdIndex,
            vrMeanIndex,
            vrStdIndex,
            meanA2Index,
            meanM2Index
        })) {
        cerr << "Invalid header in " << inputPath << endl;
        return;
    }

    map<string, Accumulator> merged;

    while (getline(fin, line)) {
        if (line.empty()) {
            continue;
        }

        vector<string> fields = splitCsvLine(line);
        int neededSize = max(meanA2Index, meanM2Index) + 1;
        if (static_cast<int>(fields.size()) < neededSize) {
            cerr << "Skipping malformed row in " << inputPath << endl;
            continue;
        }

        string vtMean = fields[vtMeanIndex];
        string vtStd = fields[vtStdIndex];
        string vrMean = fields[vrMeanIndex];
        string vrStd = fields[vrStdIndex];
        string key = keyFor(vtMean, vtStd, vrMean, vrStd);

        Accumulator& acc = merged[key];
        if (acc.count == 0) {
            acc.vtMean = vtMean;
            acc.vtStd = vtStd;
            acc.vrMean = vrMean;
            acc.vrStd = vrStd;
        }

        acc.sumA2 += stod(fields[meanA2Index]);
        acc.sumM2LogScore += stod(fields[meanM2Index]);
        acc.count++;
    }

    ofstream fout(outputPath, ios::trunc);
    fout << "vt_mean,vt_std,vr_mean,vr_std,mean_A2,mean_m2_log_score,merged_count\n";

    for (const auto& entry : merged) {
        const Accumulator& acc = entry.second;
        fout << acc.vtMean << ','
             << acc.vtStd << ','
             << acc.vrMean << ','
             << acc.vrStd << ','
             << acc.sumA2 / acc.count << ','
             << acc.sumM2LogScore / acc.count << ','
             << acc.count << '\n';
    }

    cout << inputPath << " -> " << outputPath
         << " (" << merged.size() << " merged rows)" << endl;
}

int main() {
    filesystem::create_directories("Data");

    for (int i = 1; i <= 4; i++) {
        string inputPath = "spiralSimulLogs" + to_string(i) + "/analysis_index.csv";
        string outputPath = "Data/spiralSimulLogs" + to_string(i) + "_merged.csv";
        mergeOneFile(inputPath, outputPath);
    }

    return 0;
}
