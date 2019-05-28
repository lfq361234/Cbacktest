#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include"backtest.h"
using namespace std;


int main() {
	writeAllBacktest();
	return 0;
}

////output backtest function to python in dll mode 
//#include<boost/python.hpp>
//using namespace boost::python;
//BOOST_PYTHON_MODULE(backtest)			// Python model start
//{
//	def("writeAllBacktest", writeAllBacktest);
//	def("readStrategy", readStrategy);
//	def("readParameter", readParameter);
//	def("computeCost", computeCost);
//	def("processCost", processCost);
//	def("processCRPnL", processCRPnL);
//	def("processURPnL", processURPnL);
//	def("processCPnL", processCPnL);
//	def("processMPnL", processMPnL);
//	def("writeBacktest", writeBacktest);
//	def("removeZero", removeZero);
//	def("computeMean", computeMean);
//	def("computeVariance", computeVariance);
//	def("computeSemiVariance", computeSemiVariance);
//	def("computeYR", computeYR);
//	def("computeSTD", computeSTD);
//	def("computeDownSTD", computeDownSTD);
//	def("processNAV", processNAV);
//	def("computeMDD", computeMDD);
//	def("isMarkUpDown", isMarkUpDown);
//	def("processOneWinLose", processOneWinLose);
//	def("processOneWinLoseAnalysis", processOneWinLoseAnalysis);
//	def("processAllWinLose", processAllWinLose);
//	def("writeEvaluation", writeEvaluation);
//
//}


