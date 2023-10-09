#include "pch.h"
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <locale>
#include <vector>
#include <numeric> 
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

using vec = vector<double>;
using matrix = vector<vec>;
//====================================================================================================================================================
//Option Pricing Functions
double monte_carlo_Eur_price(const int num_sims, const double S, const double K, const double r, const double v,  //Pricing European Option
							const double T, const double d, int OptionProcess, int OptionType);   
double priceExoticBySimulation(double T, double K, double S, double v, double r, double d,   //Exotic Option Pricing
							double barrier, int numberOfPaths, int process, int type);  
void calc_path_spot_prices(std::vector<double>& spot_prices, const double& r, const double& v,   //underlying underlying prices path
							const double& T, const double& d, int process);  
double payoffFunction(std::vector<double>& spot_prices, double& K, int& type, double& barrier);   //Payoff of Exotic Options
double rand_normal(double mean, double stddev);  //generate a gaussian number with different miu and sigma
double gaussian_box_muller();  //generate a gaussian number with normal distribution
//====================================================================================================================================================
//Runtime Functions
double validator(double Para, double Lbound, double Ubound);  //validation function for parameters
int validateChoice(double Choice, int numberOfChoices); // validation function for choices
matrix readCSV(string filename);  //CSV file reader, get value from the CSV file
void write(const matrix &M, int& Type, double& S, double& K, double&B, double& r,
	double& v, double& T, double& d, int& Process);  //CSV writer, return the value of Options input from CSV file
//====================================================================================================================================================
//Main function

int main()
{
	//Define exit programme indicator
	int exit = 1;
	std::cout << "Welcome to my programme!\n" << "This programme will use Monte-Carlo method to price a series of options\n\n"; 

	do
	{
		//define variables for data validation, tracking user choices and performing the calculations
		int OptionType = 0;  //OptionType definition
		double OptionValue;  //the value of the option at current time
		int OptionProcess = 0;  //OptionProcess definition 1 for BS 2 for Bachelier
		int OptionChoice = 0;  //3 kinds of modes(Single option input by hand, Single option input from CSV, Multiple options input from CSV)
		double Spot = 0;  //Spot price of the underlying
		double Strike = 0;  //strike of the option
		double r = 0;  //annualized risk-free rate
		double Vol = 0;  //annualized volatility 
		double Expiry = 0;  //the maturity of the option(year count basis)
		double d = 0;  //dividend
		double barrier = 0;  //barrier for barrier option
		int num_sims;  //Number of simulations in total

		std::cout << "Welcome to my programme!\n" << "This programme will use Monte-Carlo method to price a series of options\n";  //greatings

		//define the service user wants
		std::cout << "\nWhat kind of service do you want?";
		std::cout << "\n(1)Single Option Pricing getting data by hand\n(2)";
		std::cout << "Single Option Pricing getting data from CSV file\n(3)";
		std::cout << "Multiple Options Pricing getting data from CSV file\n";
		cin >> OptionChoice;
		validateChoice(OptionChoice, 3);

		//define the underlying price process
		std::cout << "\nPlease enter the Process of your option:";
		std::cout << "\n(1)BS(Log-normal Process)\n(2)Bachelier Process\n";
		std::cin >> OptionProcess;
		OptionProcess = validateChoice(OptionProcess, 2);

		//If user choose to pricing single option input by hand
		if (OptionChoice == 1) {

			//define the type of option and validate user inputs
			std::cout << "\nPlease enter the type of your option:";
			std::cout << "\n(1)European Call (2)European Put"; 
			std::cout << "\n(3)Asian Call (4)Asian put";
			std::cout << "\n(5)Down and In Call (6)Down and In Put";
			std::cout << "\n(7)Down and Out Call (8)Down and Out Put";
			std::cout << "\n(9)Up and In Call (10)Up and In Put";
			std::cout << "\n(11)Up and Out Call (12)Up and Out Put\n";
			std::cin >> OptionType;
			OptionType = validateChoice(OptionType, 12);

			//define Spot Price of the underlying and validate user inputs
			std::cout << "\nNow enter parameters for pricing your option:\n";
			std::cout << "Please Enter Spot Price of the underlying: $";
			std::cin >> Spot;
			Spot = validator(Spot, 0, 10000);

			//define Strike of the underlying and validate user inputs
			std::cout << "\nPlease Enter Strike of the underlying: $";
			std::cin >> Strike;
			Strike = validator(Strike, 0, 10000);

			if ((OptionType != 1) & (OptionType != 2) & (OptionType != 3) & (OptionType != 4))
			{
				std::cout << "\nPlease Enter Barrier boundary: $";
				std::cin >> barrier;
				barrier = validator(barrier, 1, 10000);
			}

			//define annualized risk free rate and validate user inputs
			std::cout << "\nPlease Enter annualized Risk-Free Rate: ";
			std::cin >> r;
			r = validator(r, -1, 1);

			//define annualized dividend rate and validate user inputs
			std::cout << "\nPlease Enter annualized Dividend Rate: ";
			std::cin >> d;
			d = validator(d, 0, 1);

			//define annualized volatility and validate user inputs
			std::cout << "\nPlease Enter annualized Volatility: ";
			std::cin >> Vol;
			Vol = validator(Vol, 0, 1);


			//define Maturity (how many years) and validate user inputs
			std::cout << "\nPlease Enter Maturity (Year Count Basis): ";
			std::cin >> Expiry;
			Expiry = validator(Expiry, 0, 1000);
		}


		if (OptionChoice == 1) 
		{
			////Calculate option price by monte carlo for vanillas
			if ((OptionType == 1) || (OptionType == 2))
			{
				num_sims = 1000000;
				OptionValue = monte_carlo_Eur_price(num_sims, Spot, Strike, r, Vol, Expiry, d, OptionProcess, OptionType);
				std::cout << "\nThe value of your European option is: " << OptionValue << std::endl;
			}
			//Calculate exotic option price with simulator function, including option type and process as input
			else
			{
				num_sims = 100000;
				OptionValue = priceExoticBySimulation(Expiry, Strike, Spot, Vol, r, d, barrier, num_sims, OptionProcess, OptionType);
				std::cout << "\nThe value of your selected Exotic Option is: " << OptionValue << std::endl;
			}
		}

		//If User choose to pricing options input from CSV file
		if ((OptionChoice == 2) || (OptionChoice == 3))
		{
			matrix M = readCSV("OptionPricing.csv");
			write(M, OptionType, Spot, Strike, barrier, r, Vol, Expiry, d, OptionProcess);
		}

		std::cout << "Do you want to price another option? Please enter 1 for yes, otherwise enter any value: ";
		std::cin >> exit;

	} while (exit == 1);  //if user press 1, run the programme again, otherwise quit

	std::cout << "Thanks for using our option pricing program. Bye!" << endl;
	return 0;
}

//====================================================================================================================================================
//validation function to check user's value parameter input and give proper value indications 
double validator(double Para, double Lbound, double Ubound)
{
	int i;
	for (i = 1; i < 5; i++)
	{
		if ((Para >= Lbound) & (Para <= Ubound) && (std::cin))
		{
			return Para;
		}
		else
		{
			//Clear console and inform user of erroneous input
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			std::cout << "You have entered wrong input, the value should be a number between " << Lbound << " and " << Ubound <<
				", please enter a new input: " << endl;
			std::cout << "(The programme will be terminated after " << 5 - i << " try)" << endl;
			std::cin >> Para;
		}
		if (i >= 5) //user can try up to 5 times
		{
			cout << "Input has failed too many times, please close program and try again." << endl;
			exit(0);
		}
	}
	return Para;

}

//Function to validate the Number of choice input and gives proper validations
int validateChoice(double Choice, int numberOfChoices)
{
	int i;
	for (i = 1; i < 5; i++)
	{
		//Check user inputs integer value within range of choices
		if ((Choice <= numberOfChoices) & (Choice > 0) && (std::cin))
		{
			return Choice;
		}
		else
		{
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			std::cout << "You have entered wrong input, you should enter the number of the choice you wish to make"
				", please enter a new input: " << endl;
			std::cout << "(The programme will be terminated after " << 5 - i << " try)" << endl;
			std::cin >> Choice;
		}
		if (i >= 5) //user can try up to 5 times
		{
			cout << "Input has failed too many times, please close program and try again." << endl;
			exit(0);
		}
	}
	return Choice;
}

//====================================================================================================================================================
//Functions to calculate the value of European Option
double monte_carlo_Eur_price(const int num_sims, const double S, const double K, const double r,
	const double v, const double T, const double d, int OptionProcess, int OptionType) {
	double S_adjust = S * exp(T * (r - d - 0.5 * v * v));
	double S_cur = 0.0;
	double payoff_sum = 0.0;
	double gauss_bm = 0.0;

	if (OptionProcess == 1) {    //BS Process(Log-normal Process or GBM process)
		for (int i = 0; i < num_sims; i++) {
			gauss_bm = gaussian_box_muller();
			S_cur = S_adjust * exp(sqrt(v * v * T) * gauss_bm);
			if (OptionType == 1) {  //European Call 
				payoff_sum += std::max(S_cur - K, 0.0);
			}
			if (OptionType == 2) {  // European Put
				payoff_sum += std::max(K - S_cur, 0.0);
			}

		}
	}

	if (OptionProcess == 2) {  //Bachelier Process
		for (int i = 0; i < num_sims; i++) {
			gauss_bm = rand_normal(0, exp(2 * (r - d) * T - 1) / (2 * (r - d)));
			S_cur = S * exp((r - d) * T) + gauss_bm * v;
			if (OptionType == 1) {  //European Call 
				payoff_sum += std::max(S_cur - K, 0.0);
			}
			if (OptionType == 2) {  //European Put 
				payoff_sum += std::max(K - S_cur, 0.0);
			}
		}
	}

	return (payoff_sum / static_cast<double>(num_sims)) * exp(-r * T);
}

//Functions to calculate the value of Exotic Options
double priceExoticBySimulation(double T, double K, double S, double v, double r, double d,
	double barrier, int numberOfPaths, int process, int type)
{
	double payoff_sum;
	payoff_sum = 0;
	std::vector<double> spot_prices(250); // Initiate spot price vector and assign initial spot value
	spot_prices[0] = S;
	for (int i = 0; i < numberOfPaths; i++) {
		calc_path_spot_prices(spot_prices, r, v, T, d, process); //Generate spot prices along path 
		payoff_sum += payoffFunction(spot_prices, K, type, barrier); //Calculate payoff for specific exotic option
	}
	double discount_payoff_avg = (payoff_sum / numberOfPaths) * exp(-r * T);
	return discount_payoff_avg;

}

// Vector of spot prices to be filled in
void calc_path_spot_prices(std::vector<double> & spot_prices,  // Vector of spot prices to be filled in
	const double& r,   // Risk free interest rate (constant)
	const double& v,   // Volatility of underlying (constant)
	const double& T,  // Expiry
	const double& d,
	int process) { // dividend
// Since the drift and volatility of the asset are constant
// we will precalculate as much as possible for maximum efficiency
	double dt = T / static_cast<double>(spot_prices.size());
	double drift = exp(dt * (r - d - 0.5 * v * v));
	double vol = sqrt(v * v * dt);
	
	//Fill spot price vector with generated values for selected process
	for (size_t i = 1; i < spot_prices.size(); i++) {
		if (process == 1) {
			double gauss_bm = gaussian_box_muller();
			spot_prices[i] = spot_prices[i - 1] * drift * exp(vol * gauss_bm);
		}

		if (process == 2) {
			double gauss_bm = rand_normal(0, exp(2 * (r - d) * (dt * i) - 1) / (2 * (r - d)));
			spot_prices[i] = spot_prices[0] * exp((r - d) * (dt * i)) + gauss_bm * v;
		}
	}
}

//Payoff of Exotic Options
double payoffFunction(std::vector<double> & spot_prices, double& K, int& type, double& barrier) {
	double payoff;
	double price = 0.0;
	double spot = 0.0;
	int num_times = spot_prices.size();
	double sum = 0;
	double arith_mean = 0;
	double in = 0;
	double out = 0;
	if (type == 3) //Calculate payoff for Asian Call
	{
		sum = std::accumulate(spot_prices.begin(), spot_prices.end(), 0);
		arith_mean = sum / num_times;
		payoff = arith_mean - K;
		price = std::max(payoff, 0.0);
	}
	if (type == 4) //Calculate payoff for Asian Put
	{
		sum = std::accumulate(spot_prices.begin(), spot_prices.end(), 0);
		arith_mean = sum / num_times;
		payoff = K - arith_mean;
		price = std::max(payoff, 0.0);
	}
	if (type == 5) //Calculate payoff for Down and In Call
	{
		for (size_t i = 0; i < spot_prices.size(); i++)
		{
			spot = spot_prices[i];
			if (spot < barrier)
			{
				in = 1;
			}
		}
		if (in == 1)
		{
			payoff = spot_prices[num_times - 1] - K;
			price = std::max(payoff, 0.0);
		}
	}
	if (type == 6) //Calculate payoff for Down and In Put
	{
		for (size_t i = 0; i < spot_prices.size(); i++)
		{
			if (spot_prices[i] < barrier)
			{
				in = 1;
			}
		}
		if (in == 1)
		{
			payoff = K - spot_prices[num_times - 1];
			price = std::max(payoff, 0.0);
		}
	}
	if (type == 7) //Calculate payoff for Down and Out Call
	{
		for (size_t i = 0; i < spot_prices.size(); i++)
		{
			if (spot_prices[i] < barrier)
			{
				out = 1;
			}
		}

		if (out != 1)
		{
			payoff = spot_prices[num_times - 1] - K;
			price = std::max(payoff, 0.0);
		}
	}
	if (type == 8)//Calculate payoff for Down and Out Put
	{
		for (size_t i = 0; i < spot_prices.size(); i++)
		{
			if (spot_prices[i] < barrier)
			{
				out = 1;
			}
		}
		if (out != 1)
		{
			payoff = K - spot_prices[num_times - 1];
			price = std::max(payoff, 0.0);
		}
	}
	if (type == 9) //Calculate payoff for Up and In Call
	{
		for (size_t i = 0; i < spot_prices.size(); i++)
		{
			if (spot_prices[i] > barrier)
			{
				in = 1;
			}
		}
		if (in == 1)
		{
			payoff = spot_prices[num_times - 1] - K;
			price = std::max(payoff, 0.0);
		}
	}
	if (type == 10) //Calculate payoff for Up and In Put
	{
		for (size_t i = 0; i < spot_prices.size(); i++)
		{
			if (spot_prices[i] > barrier)
			{
				in = 1;
			}
		}
		if (in == 1)
		{
			payoff = K - spot_prices[num_times - 1];
			price = std::max(payoff, 0.0);
		}
	}
	if (type == 11) //Calculate payoff for Up and Out Call
	{
		for (size_t i = 0; i < spot_prices.size(); i++)
		{
			if (spot_prices[i] > barrier)
			{
				out = 1;
			}
		}

		if (out != 1)
		{
			payoff = spot_prices[num_times - 1] - K;
			price = std::max(payoff, 0.0);
		}

	}
	if (type == 12) //Calculate payoff for Up and Out Put
	{
		for (size_t i = 0; i < spot_prices.size(); i++)
		{
			if (spot_prices[i] > barrier)
			{
				out = 1;
			}
		}

		if (out != 1)
		{
			payoff = K - spot_prices[num_times - 1];
			price = std::max(payoff, 0.0);
		}

	}

	return price;
}

//====================================================================================================================================================
// A simple implementation of the Box-Muller algorithm, used to generate standard normal gaussian random numbers
double gaussian_box_muller() {
	double x = 0.0;
	double y = 0.0;
	double euclid_sq = 0.0;

	// Continue generating two uniform random variables
	// until the square of their "euclidean distance" 
	// is less than unity
	do {
		x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		euclid_sq = x * x + y * y;
	} while (euclid_sq >= 1.0);

	return x * sqrt(-2 * log(euclid_sq) / euclid_sq);
}

//generate a normal gaussian number with different miu and sigma
double rand_normal(double mean, double stddev) 
{//Box muller method
	static double n2 = 0.0;
	static int n2_cached = 0;
	if (!n2_cached)
	{
		double x, y, r;
		do
		{
			x = 2.0*rand() / RAND_MAX - 1;
			y = 2.0*rand() / RAND_MAX - 1;

			r = x * x + y * y;
		} while (r == 0.0 || r > 1.0);
		{
			double d = sqrt(-2.0*log(r) / r);
			double n1 = x * d;
			n2 = y * d;
			double result = n1 * stddev + mean;
			n2_cached = 1;
			return result;
		}
	}
	else
	{
		n2_cached = 0;
		return n2 * stddev + mean;
	}
}

//==================================================================================================================================================================
//CSV reader, get the value from CSV file
matrix readCSV(string filename)
{
	matrix M;

	ifstream in(filename);  //to check whether the file is open or not
	if (!in.is_open()) {
		std::exit(EXIT_FAILURE);
	}

	//ignore the first 2 lines of notes
	//ignore the first 500 characters, or until first \n, whichever is met first so that we can ignore the first line of the CSV file
	in.ignore(500, '\n'); 
	in.ignore(500, '\n');

	string line;
	while (getline(in, line))                   // read a whole line of the file
	{
		stringstream ss(line);                     // put it in a stringstream (internal stream)
		vec row;
		string data;
		while (getline(ss, data, ','))           // read (string) items up to a comma
		{
			row.push_back(stod(data));            // use stod() to convert to double; put in row vector
		}
		if (row.size() > 0) M.push_back(row);    // add non-empty rows to matrix
	}

	return M;
}

//CSV writer, return the value of Options input from CSV file
void write(const matrix &M, int& Type, double& S, double& K, double&B, double& r,
 double& v, double& T, double& d, int& Process)
{
	double portfolioPrice = 0;
	double portfolioPrice_ = 0;
	// cout << fixed;
	for (auto row : M)
	{
		for (int i = 0; i < row.size() ; i = i + 8)
		{
			Type = row[i];
			S = row[i+1];
			K = row[i+2];
			B = row[i+3];
			r = row[i+4];
			d = row[i+5];
			v = row[i+6];
			T = row[i+7];

			if ((Type == 1) || (Type == 2))
			{
				int num_sims = 1000000;
				portfolioPrice_ = monte_carlo_Eur_price(num_sims, S, K, r, v, T, d, Process, Type);
			}

			//Calculate exotic option price with simulator function, including option type and process as input
			if ((Type != 1) & (Type != 2))
			{
				int num_sims = 100000;
				portfolioPrice_ = priceExoticBySimulation(T, K, S, v, r, d, B, num_sims, Process, Type);
			}

			portfolioPrice += portfolioPrice_;
		}
	}

	std::cout << "\nThe value of your portfolio is: " << portfolioPrice << std::endl;
}
