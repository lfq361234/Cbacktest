#pragma once
#include<vector>
#include<iostream>
//lookup backtest.cpp for illustration
std::vector<std::pair<std::string,std::vector<std::pair<double, double>>>>  readStrategy();
void readParameter(double&, double&, double&, double&, double&);
double computeCost(const double, const double, const double, const double);
std::vector<double> processCost(const std::vector<std::pair<std::string, std::vector<std::pair<double, double>>>>&, const double, const double);
std::vector<double> processCRPnL(const std::vector<std::pair<std::string, std::vector<std::pair<double, double>>>>&, const std::vector<double>&);
std::vector<double> processURPnL(const std::vector<std::pair<std::string, std::vector<std::pair<double, double>>>>&, const double, const double);
std::vector<double> processCPnL(const std::vector<double>&, const std::vector<double>&);
std::vector<double> processMPnL(const std::vector<double>&);
void writeBacktest(const std::vector<std::pair<std::string, std::vector<std::pair<double, double>>>>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);
std::vector<double> removeZero(const std::vector<double>&);
double computeMean(const std::vector<double>&);
double computeVariance(const std::vector<double>&);
double computeSemiVariance(const std::vector<double>&);
double computeYR(const std::vector<double>&, const double, const double);
double computeSTD(const std::vector<double>&, const double, const double);
double computeDownSTD(const std::vector<double>&, const double, const double, const double);
std::vector<double> processNAV(const std::vector<double>&, const double);
double computeMDD(const std::vector<double>&);
bool isMarkUpDown(const std::vector<double>&);
std::vector<std::pair<std::string, double>> processOneWinLose(const int, const std::vector<std::pair<std::string, std::vector<std::pair<double, double>>>>&, const double, const double);
std::vector<double> processOneWinLoseAnalysis(const int, const std::vector<std::pair<std::string, std::vector<std::pair<double, double>>>>&, const double, const double);
std::vector<std::vector<std::pair<std::string, double>>> processAllWinLose(const std::vector<std::pair<std::string, std::vector<std::pair<double, double>>>>&, const double, const double);
void writeEvaluation(const std::vector<std::pair<std::string, std::vector<std::pair<double, double>>>>&, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double, const std::vector<std::vector<std::pair<std::string, double>>>&);
void writeAllBacktest();