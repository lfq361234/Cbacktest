#include"backtest.h"
#include<fstream>
#include<sstream>
#include<numeric>
#include<limits>
using namespace std;
/**@brief: read strategy csv data from "form_strategy.csv"
	  csv data has to follow the format
	  first row:time,asset1,strategy1,asset2,strategy2,...
	  other rows: dataPoint
	  number of asset and strategy should be same

	  @param_out: use vector to store data because it is fast for access
	  vector of pairs<time,assetStrategyPairs>
	  assetStrategyPairs: vector of pairs<assetprice, strategy>
	  if file can not open or the number of columns is wrong, the output will be empty

	 examples:
	  initial time, first asset, price
	  strategyData[0].second[0].first
	  initial time, first asset, strategy
	  strategyData[0].second[0].second
	  initial time, second asset, strategy
	  strategyData[0].second[1].second
**/
vector<pair<string, vector<pair<double, double>>>>  readStrategy() {
	vector<pair<string, vector<pair<double, double>>>> strategyData;
	fstream file;
	file.open("./input_file/form_strategy.csv");
	if (file.is_open()) {
		//first row does not contain content of data, we use it to check proper number of columns
		string line;
		getline(file, line);
		int columnCount = 0;
		istringstream firstRow(line);
		string s;
		while (getline(firstRow, s, ',')) {
			columnCount++;
		}
		// if number of columns is right, start to read data
		//output format: vector of pairs<time,assetStrategyPairs>
		//assetStrategyPairs: vector of pairs<asset, strategy>
		if (columnCount % 2 == 1) {
			while (getline(file, line)) {
				istringstream tempLine(line);
				//read time data
				string dataTime;
				getline(tempLine, dataTime, ',');
				//read other data and store in assetStrategyPairs
				vector<pair<double, double>> assetStrategyPairs;
				string s1, s2;
				while (getline(tempLine, s1, ',') && getline(tempLine, s2, ',')) {
					double d1 = atof(s1.c_str());
					double d2 = atof(s2.c_str());
					pair<double, double> p1 = make_pair(d1, d2);
					assetStrategyPairs.push_back(p1);
				}
				strategyData.push_back(make_pair(dataTime, assetStrategyPairs));
			}
		}
		else {
			cout << "number of columns is wrong, so return vector from parser is empty" << endl;
		}
	}
	else {
		cout << "can not open csv file, so return vector from parser is empty" << endl;
	}
	file.close();
	return strategyData;
}


//read "parameter.csv" for setting parameters
//first row: eachTradeCost,percetageOfTradeVolumnCost,dataFreq,riskFreeRate,investMoney
//second row: data
void readParameter(double& eachTradeCost, double& percetageOfTradeVolumnCost, double& dataFreq, double& riskFreeRate, double& investMoney) {
	fstream file;
	file.open("./input_file/parameter.csv");
	if (file.is_open()) {
		//first row dose not contain data
		string line;
		getline(file, line);
		//read data from second row 
		string data;
		//first parameter
		getline(file, data, ',');
		eachTradeCost = atof(data.c_str());
		//second parameter ...
		getline(file, data, ',');
		percetageOfTradeVolumnCost = atof(data.c_str());
		getline(file, data, ',');
		dataFreq = atof(data.c_str());
		getline(file, data, ',');
		riskFreeRate = atof(data.c_str());
		getline(file, data, ',');
		investMoney = atof(data.c_str());

	}
	else {
		cout << "can not open parameter.csv" << endl;
	}
	file.close();
}

//basic cost formula
double computeCost(const double price, const double strategy, const double eachTradeCost, const double percetageOfTradeVolumnCost) {
	if (strategy == 0) {
		return 0;
	}
	else {
		return (eachTradeCost + abs(strategy) * price * percetageOfTradeVolumnCost) * -1;
	}
}

//calculate total cost of all assets at each time point
vector<double> processCost(const vector<pair<string, vector<pair<double, double>>>>& strategyData, const double eachTradeCost, const double percetageOfTradeVolumnCost) {
	vector<double> result;
	for (int i = 0; i < strategyData.size(); i++) {
		double momentCost = 0;
		for (int j = 0; j < strategyData[i].second.size(); j++) {
			momentCost = momentCost + computeCost(strategyData[i].second[j].first, strategyData[i].second[j].second, eachTradeCost, percetageOfTradeVolumnCost);
		}
		result.push_back(momentCost);
	}
	return result;
}

//compute cumulative realized PnL with transaction cost for all assets during whole time period
vector<double> processCRPnL(const vector<pair<string, vector<pair<double, double>>>>& strategyData, const vector<double>& transactionCost) {
	double CTcost = 0;
	double CPnL = 0;
	double Csum = 0;
	vector<double> result;
	for (int i = 0; i < strategyData.size(); i++) {
		//cumulative cost
		CTcost = CTcost + transactionCost[i];
		//cumulative PnL without cost
		double momentPnL = 0;
		for (int j = 0; j < strategyData[i].second.size(); j++) {
			momentPnL = momentPnL + strategyData[i].second[j].first * strategyData[i].second[j].second;
		}
		CPnL = CPnL + momentPnL;
		//cumulative PnL with cost
		Csum = CTcost + CPnL;
		result.push_back(Csum);
	}
	return result;
}

// compute unrealized PnL with transaction cost for all assets  at each moment
vector<double> processURPnL(const vector<pair<string, vector<pair<double, double>>>>& strategyData, const double eachTradeCost, const double percetageOfTradeVolumnCost) {
	// initialization of sum of strategy per asset 
	vector<double> sumStrategyPerAsset;
	if (!strategyData.empty()) {
		for (int j = 0; j < strategyData[0].second.size(); j++) {
			sumStrategyPerAsset.push_back(0);
		}
	}
	vector<double> result;
	for (int i = 0; i < strategyData.size(); i++) {
		double momentURPnL = 0;
		//calcuate sum of strategy per asset
		for (int j = 0; j < strategyData[i].second.size(); j++) {
			sumStrategyPerAsset[j] += strategyData[i].second[j].second;
			//add moment unrealized PnL per asset to momentURPnL 
			if (sumStrategyPerAsset[j] != 0) {
				momentURPnL += strategyData[i].second[j].first * -1 * sumStrategyPerAsset[j] + computeCost(strategyData[i].second[j].first, -1 * sumStrategyPerAsset[j], eachTradeCost, percetageOfTradeVolumnCost);
			}
		}
		result.push_back(momentURPnL);
	}
	return result;
}

//compute cumulative PnL with realized PnL and unrealized PnL
vector<double> processCPnL(const vector<double>& CRPnL, const vector<double>& URPnL) {
	vector<double> result;
	if (CRPnL.size() == URPnL.size()) {
		for (int i = 0; i < CRPnL.size(); i++) {
			result.push_back(CRPnL[i] + URPnL[i]);
		}
	}
	else {
		cout << "wrong input in calculation of process CPnL" << endl;
	}
	return result;
}

//compute momentPnL with cumulative PnL
vector<double>processMPnL(const vector<double>& CPnL) {
	vector<double>result;
	for (int i = 0; i < CPnL.size(); i++) {
		if (i == 0) {
			result.push_back(CPnL[i]);
		}
		else {
			result.push_back(CPnL[i] - CPnL[i - 1]);
		}
	}
	return result;
}

//output format: time,transactionCost, realizedCPnL, unrealizedPnL, cumulativePnL, momentPnL
void writeBacktest(const vector<pair<string, vector<pair<double, double>>>>& strategyData, const vector<double>& TCost, const vector<double>& CRPnL, const vector<double>& URPnL, const vector<double>& CPnL, const vector<double>& MPnL) {
	remove("./output_file/backtest.csv");
	fstream file;
	file.open("./output_file/backtest.csv", ios::out);
	if (file.is_open()) {
		file << "time,transactionCost,realizedCPnL,unrealizedPnL,cumulativePnL,momentPnL" << endl;
		for (int i = 0; i < strategyData.size(); i++) {
			file << strategyData[i].first << "," << TCost[i] << "," << CRPnL[i] << "," << URPnL[i] << "," << CPnL[i] << "," << MPnL[i] << endl;
		}
	}
	else {
		cout << "can not write into backtest.csv" << endl;
	}
	file.close();
}

//eliminate 0 from momentPnL for holdingPnL
//eliminate 0 from strategy for judging markup or markdown
vector<double> removeZero(const vector<double>& v) {
	vector<double> result;
	for (int i = 0; i < v.size(); i++) {
		if (v[i] != 0) {
			result.push_back(v[i]);
		}
	}
	return result;
}

double computeMean(const vector<double>& v) {
	if (v.empty()) return numeric_limits<double>::quiet_NaN();
	return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double computeVariance(const vector<double>& v) {
	if (v.size() <= 1) return numeric_limits<double>::quiet_NaN();
	double sum = 0;
	double sumSqu = 0;
	for (int i = 0; i < v.size(); i++) {
		sum = sum + v[i];
		sumSqu = sumSqu + pow(v[i], 2);
	}
	return sumSqu / (v.size() - 1) - v.size() * pow(sum / v.size(), 2) / (v.size() - 1);
}

double computeSemiVariance(const vector<double>& v) {
	if (v.size() <= 1) return numeric_limits<double>::quiet_NaN();
	double mean = computeMean(v);
	double numerator = 0;
	for (int i = 0; i < v.size(); i++) {
		if (v[i] < mean) {
			numerator += pow(v[i] - mean, 2);
		}
	}
	return numerator / (v.size() - 1);
}

double computeYR(const vector<double>& holdingPnL, const double dataFreq, const double investMoney) {
	if (computeMean(holdingPnL) == numeric_limits<double>::quiet_NaN()) return numeric_limits<double>::quiet_NaN();
	return computeMean(holdingPnL) * dataFreq / investMoney;
}

double computeSTD(const vector<double>& holdingPnL, const double dataFreq, const double investMoney) {
	if (computeVariance(holdingPnL) == numeric_limits<double>::quiet_NaN()) return numeric_limits<double>::quiet_NaN();
	return sqrt(computeVariance(holdingPnL)) * sqrt(dataFreq) / investMoney;
}

//Use minimum acceptable return to calculate target average daily PnL,and downside STD
//Assume minimum acceptable return is equal to risk-free rate in parameter setting
double computeDownSTD(const vector<double>& holdingPnL, const double dataFreq, const double investMoney, const double MAC) {
	if (holdingPnL.size() <= 1) return numeric_limits<double>::quiet_NaN();
	double targetMean = MAC * investMoney / dataFreq;
	double numerator = 0;
	for (int i = 0; i < holdingPnL.size(); i++) {
		if (holdingPnL[i] < targetMean) {
			numerator += pow(holdingPnL[i] - targetMean, 2);
		}
	}
	return sqrt(numerator / holdingPnL.size()) * sqrt(dataFreq) / investMoney;
}
//get netAssetValue form investMoney and cumulativePnL
vector<double> processNAV(const vector<double>& cumulativePnL, const double investMoney) {
	vector<double> result;
	for (int i = 0; i < cumulativePnL.size(); i++) {
		result.push_back(investMoney + cumulativePnL[i]);
	}
	return result;
}

//compute MaxDrawDown from Net Asset Value
double computeMDD(const vector<double>& NAV) {
	double peak = -99999;
	double DD = 0;
	double MDD = 0;
	for (int i = 0; i < NAV.size(); i++) {
		if (NAV[i] > peak) {
			peak = NAV[i];
		}
		DD = (peak - NAV[i]) / peak;
		if (DD > MDD) {
			MDD = DD;
		}
	}
	return MDD;
}

// Judging if there is markup or markdown from single asset strategy after removing zeros 
bool isMarkUpDown(const vector<double>& strategy) {
	double countFalsePairs = 0;
	if (strategy.empty()) return false;
	else {
		for (int i = 0; i < strategy.size(); i = i + 2) {
			if (i + 1 < strategy.size()) {
				if (strategy[i] + strategy[i + 1] == 0) {
					countFalsePairs++;
				}
			}
		}
	}
	if (countFalsePairs == static_cast<double>(strategy.size()) / 2 || countFalsePairs == static_cast<double>(strategy.size() - 1) / 2) {
		return false;
	}
	else {
		return true;
	}
}
//deal with win/lose  to get vector<pair<string, double>> 
//vector of pairs(transaction time, win/lose amount)
//only for one asset without markup/markdown
vector<pair<string, double>> processOneWinLose(const int assetIndex, const vector<pair<string, vector<pair<double, double>>>>& strategyData, const double eachTradeCost, const double percetageOfTradeVolumnCost) {
	//get the strategy data for assetIndex
	vector<double> strategy;
	if (!strategyData.empty()) {
		for (int i = 0; i < strategyData.size(); i++) {
			if (assetIndex < strategyData[i].second.size()) {
				strategy.push_back(strategyData[i].second[assetIndex].second);
			}
		}
	}
	// find out the moment of transaction for assetIndex limited to no markup/markdown
	vector<int> indexRecord;
	if (!isMarkUpDown(removeZero(strategy))) {
		for (int i = 0; i < strategyData.size(); i++) {
			// if strategy of assetIndex is not equal to 0, record it to indexRecord
			if (assetIndex < strategyData[i].second.size()) {
				if (strategyData[i].second[assetIndex].second != 0) {
					indexRecord.push_back(i);
				}
			}
		}
	}
	//fillout winLoseTime, winLoseAmount for assetIndex
	vector<pair<string, double>> winLose;
	for (int j = 0; j < indexRecord.size(); j = j + 2) {
		if (j + 1 < indexRecord.size()) {
			if (indexRecord[j + 1] < strategyData.size() && indexRecord[j] < strategyData.size()) {
				if (assetIndex < strategyData[indexRecord[j + 1]].second.size() && assetIndex < strategyData[indexRecord[j]].second.size()) {
					string winLoseTime = strategyData[indexRecord[j + 1]].first;
					double p2 = strategyData[indexRecord[j + 1]].second[assetIndex].first;
					double q2 = strategyData[indexRecord[j + 1]].second[assetIndex].second;
					double c2 = computeCost(p2, q2, eachTradeCost, percetageOfTradeVolumnCost);
					double p1 = strategyData[indexRecord[j]].second[assetIndex].first;
					double q1 = strategyData[indexRecord[j]].second[assetIndex].second;
					double c1 = computeCost(p1, q1, eachTradeCost, percetageOfTradeVolumnCost);
					double winLoseAmount = p2 * q2 + c2 + p1 * q1 + c1;
					winLose.push_back(make_pair(winLoseTime, winLoseAmount));
				}
			}
		}
	}
	return winLose;
}

//process win/lose analysis for one asset without markup/markdown
//0,1,2,3,4,5
//no. of transaction,no. of win,total amount of win,no. of lose,total amount of lose,win rate
vector<double> processOneWinLoseAnalysis(const int assetIndex, const vector<pair<string, vector<pair<double, double>>>>& strategyData, const double eachTradeCost, const double percetageOfTradeVolumnCost) {
	vector<pair<string, double>> winLose= processOneWinLose(assetIndex,strategyData,eachTradeCost,percetageOfTradeVolumnCost);
	vector<double> result;
	if (!winLose.empty()) {
		double wincount = 0;
		double losecount = 0;
		double winAmount = 0;
		double loseAmount = 0;
		for (int i = 0; i < winLose.size(); i++) {
			if (winLose[i].second > 0) {
				wincount++;
				winAmount += winLose[i].second;
			}
			else {
				losecount++;
				loseAmount+= winLose[i].second;
			}
		}
		result.push_back(wincount + losecount);
		result.push_back(wincount);
		result.push_back(winAmount);
		result.push_back(losecount);
		result.push_back(loseAmount);
		result.push_back(wincount / (wincount + losecount));
	}
	return result;
}

//return the report of win/lose analysis: allWinLose
//allWinLose[i][j].first: the time moment of the ith asset in strategtData and the jth winlose of the strategy
//allWinLose[i][j].second: the win/lose amount of the ith asset in strategtData and the jth winlose of the strategy
vector<vector<pair<string, double>>> processAllWinLose(const vector<pair<string, vector<pair<double, double>>>>& strategyData, const double eachTradeCost, const double percetageOfTradeVolumnCost) {
	vector<vector<pair<string, double>>> result;
	if (!strategyData.empty()) {
		//A represents index of each asset
		for (int A = 0; A < strategyData[0].second.size(); A++) {
			//get the strategy of each asset
			vector<double> strategy;
			for (int i = 0; i < strategyData.size(); i++) {
				if (A < strategyData[i].second.size()) {
					strategy.push_back(strategyData[i].second[A].second);
				}
			}
			// If strategy does not contain markup or markdown
			if (!isMarkUpDown(removeZero(strategy))) {
				vector<pair<string, double>> winLosePerAsset = processOneWinLose(A, strategyData, eachTradeCost, percetageOfTradeVolumnCost);
				result.push_back(winLosePerAsset);
			}
			// If strategy  contains markup or markdown, push nothing in result
			else {
				vector<pair<string, double>> winLosePerAsset;
				result.push_back(winLosePerAsset);
			}
		}
	}
	return result;
}
void writeEvaluation(const vector<pair<string, vector<pair<double, double>>>>& strategyData, const double eachTradeCost, const double percetageOfTradeVolumnCost, const double dataFreq, const double riskFreeRate, const double investMoney, const double totalPnL, const double annualReturn, const double annualSTD, const double annualDownSTD, const double sharpeRatio, const double sortinoRatio, const double maxDrawDown, const vector<vector<pair<string, double>>>& allWinLose) {
	remove("./output_file/evaluation.csv");
	fstream file;
	file.open("./output_file/evaluation.csv", ios::out);
	if (file.is_open()) {
		file << "start time,end time,eachTradeCost,percetageOfTradeVolumnCost,dataFreq,riskFreeRate,investMoney" << endl;
		if (strategyData.size() >= 2) {
			file << strategyData[0].first << "," << strategyData[strategyData.size() - 1].first << "," << eachTradeCost << "," << percetageOfTradeVolumnCost << "," << dataFreq << "," << riskFreeRate << "," << investMoney << endl;
		}
		file << endl;
		file << "totalPnL,annualReturn,annualSTD,annualDownSTD,Sharpe Ratio,Sortino Ratio,Max DrawDown" << endl;
		file << totalPnL << "," << annualReturn << "," << annualSTD << "," << annualDownSTD << "," << sharpeRatio << "," << sortinoRatio << "," << maxDrawDown << endl;
		file << endl;
		for (int i = 0; i < allWinLose.size(); i++) {
			file << "Asset Order," << i << endl;
			file << "Time,Realized Win/Lose, Amount" << endl;
			if (allWinLose[i].empty()) {
				file << "unanalyzable" << endl;
			}
			else {
				double winNum = 0;
				double loseNum = 0;
				double winAmount = 0;
				double loseAmount = 0;
				for (int j = 0; j < allWinLose[i].size(); j++) {
					if (allWinLose[i][j].second > 0) {
						file << allWinLose[i][j].first << ",Win," << allWinLose[i][j].second << endl;
						winNum++;
						winAmount += allWinLose[i][j].second;
					}
					else {
						file << allWinLose[i][j].first << ",Lose," << allWinLose[i][j].second << endl;
						loseNum++;
						loseAmount += allWinLose[i][j].second;
					}
				}
				file << "no. of transaction,no. of win,total amount of win,no. of lose,total amount of lose,win rate" << endl;
				file << winNum + loseNum << "," << winNum << "," << winAmount << "," << loseNum << "," << loseAmount << "," << winNum / (winNum + loseNum) << endl;
			}
			file << endl;
		}
	}
	else {
		cout << "can not write into evaluation.csv" << endl;
	}
	file.close();
}
//conduct backtesting for trading strategy
//input: "form_strategy.csv", "parameter.csv"
//output: "backtest.csv","evaluation.csv"
//readme describe format of input and output
void writeAllBacktest() {
	//read data from "form_strategy.csv"
	vector<pair<string, vector<pair<double, double>>>> strategyData = readStrategy();
	//read parameter setting in "parameter.csv"
	double eachTradeCost = numeric_limits<double>::quiet_NaN();
	double percetageOfTradeVolumnCost = numeric_limits<double>::quiet_NaN();
	double dataFreq = numeric_limits<double>::quiet_NaN();
	double riskFreeRate = numeric_limits<double>::quiet_NaN();
	double investMoney = numeric_limits<double>::quiet_NaN();
	readParameter(eachTradeCost, percetageOfTradeVolumnCost, dataFreq, riskFreeRate, investMoney);
	//calculte cost and PnL data
	vector<double>transactionCost = processCost(strategyData, eachTradeCost, percetageOfTradeVolumnCost);
	vector<double>realizedCPnL = processCRPnL(strategyData, transactionCost);
	vector<double>unrealizedPnL = processURPnL(strategyData, eachTradeCost, percetageOfTradeVolumnCost);
	vector<double>cumulativePnL = processCPnL(realizedCPnL, unrealizedPnL);
	vector<double>momentPnL = processMPnL(cumulativePnL);
	//output to backtest.csv
	writeBacktest(strategyData, transactionCost, realizedCPnL, unrealizedPnL, cumulativePnL, momentPnL);

	//calculate evaluation data based on PnL 
	double totalPnL = numeric_limits<double>::quiet_NaN();
	if (!cumulativePnL.empty()) {
		totalPnL = cumulativePnL[cumulativePnL.size() - 1];
	}
	vector<double> holdingPnL = removeZero(momentPnL);
	double annualReturn = computeYR(holdingPnL, dataFreq, investMoney);
	double annualSTD = computeSTD(holdingPnL, dataFreq, investMoney);
	double annualDownSTD = computeDownSTD(holdingPnL, dataFreq, investMoney, riskFreeRate);
	double sharpeRatio = numeric_limits<double>::quiet_NaN();
	double sortinoRatio = numeric_limits<double>::quiet_NaN();
	if (annualReturn != numeric_limits<double>::quiet_NaN() && riskFreeRate != numeric_limits<double>::quiet_NaN() && annualSTD != numeric_limits<double>::quiet_NaN() && annualDownSTD != numeric_limits<double>::quiet_NaN() && annualSTD != 0 && annualDownSTD != 0) {
		sharpeRatio = (annualReturn - riskFreeRate) / annualSTD;
		sortinoRatio = (annualReturn - riskFreeRate) / annualDownSTD;
	}
	vector<double> netAssetValue = processNAV(cumulativePnL, investMoney);
	double maxDrawDown = computeMDD(netAssetValue);
	vector<vector<pair<string, double>>> allWinLose = processAllWinLose(strategyData, eachTradeCost, percetageOfTradeVolumnCost);
	//output to evaluation.csv
	writeEvaluation(strategyData, eachTradeCost, percetageOfTradeVolumnCost, dataFreq, riskFreeRate, investMoney, totalPnL, annualReturn, annualSTD, annualDownSTD, sharpeRatio, sortinoRatio, maxDrawDown, allWinLose);

}