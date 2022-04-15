#ifndef FLUEDED_H
#define FLUEDED_H
#define _USE_MATH_DEFINES
#endif // FLUEDED_H

#pragma comment(linker, "/STACK:128000000")

#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <iomanip>
#include <iterator>
#include <pqxx/pqxx>
#include <cstdio>
#include <windows.h>
#include <thread>
#include <mutex>
#include "Pora.h"
#include "supportFunc.h"
#include "pSql.h"
#include "SimpleTimer.h"

std::mutex mtx;
//using namespace std;
//int summmm = 0;

void Generation_Poras_gauss(double Rad, double core, std::vector<Pora>& poras, double N, double sii) {
	std::random_device rd_R;
	std::normal_distribution<> uid_R(Rad, sii);
	bool kj = 0;
	double Rad_tmp = 0;
	double pot_por_rad = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     uid_R(rd_R)
				poras.push_back(Pora(Rad * core + 1.0 * i * (Rad + Rad) * core, Rad * core + 1.0 * j * (Rad + Rad) * core, Rad * core + 1.0 * k * (Rad + Rad) * core, uid_R(rd_R), 0, 0, 0, 0, 0));
			}
		}
	}
}
void Generation_Poras_experiment(std::vector<double>& x, std::vector<double>& y, std::vector<double>& interp_k, std::vector<double>& interp_b, double Rad, double core, std::vector<Pora>& poras, double N) {
	double min_x = *min_element(x.begin(), x.end());
	double max_x = *max_element(x.begin(), x.end());
	double min_y = *min_element(y.begin(), y.end());
	double max_y = *max_element(y.begin(), y.end());
	std::random_device rd_X;
	std::uniform_real_distribution<> uid_X(min_x + 0.1, max_x - 0.1);
	std::random_device rd_Y;
	std::uniform_real_distribution<> uid_Y(min_y + 0.1, max_y - 0.1);
	bool kj = 0;
	double Rad_tmp = 0;
	double pot_por_rad = 0;
	int start_time = clock();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				kj = 0;
				while (!kj) {
					pot_por_rad = uid_X(rd_X);
					if (uid_Y(rd_Y) <= linear(x, y, interp_k, interp_b, pot_por_rad, x.size())) {
						Rad_tmp = pot_por_rad;
						kj = 1;
					}
				}
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     Rad_tmp
				poras.push_back(Pora(Rad * core + 1.0 * i * (Rad + Rad) * core, Rad * core + 1.0 * j * (Rad + Rad) * core, Rad * core + 1.0 * k * (Rad + Rad) * core, Rad_tmp, 0, 0, 0, 0, 0));
			}
		}
	}
}

int rabotas_in(std::vector<Pora>& poras, int i, double pres, double DSIGMA, double SIGMA) {
	double SM = 0;
	double SMZ = 0;
	double surf_empty = 0;
	double nei_zap = 0;
	int qq = 0;
	double work;
	for (int j = 0; j < poras[i].get_V_neighbor_Num().size(); j++) {

		SM += +poras[i].get_neighbor_Area(j);
		if (poras[poras[i].get_neighbor_Num(j)].get_filled()) {
			nei_zap++;
			SMZ += poras[i].get_neighbor_Area(j);
		}
		if (!(poras[poras[i].get_neighbor_Num(j)].get_filled())) {
			surf_empty += poras[i].get_neighbor_Area(j);
		}
	}
	work = -pres * 4.0 * pow(10.0, -9.0) * M_PI / 3.0 * pow(poras[i].get_Rad(), 3.0) +
		(4.0 * M_PI * poras[i].get_Rad() * poras[i].get_Rad() - SM) * DSIGMA +
		SIGMA * (surf_empty - SMZ);
	return (work < 0);
}

double rabotas_out(std::vector<Pora>& poras, int number, double pres, double DSIGMA, double SIGMA) {
	double SM = 0;
	double SMZ = 0;
	double surf_empty = 0;
	double nei_zap = 0;
	double work = 0;
	for (int j = 0; j < poras[number].get_V_neighbor_Num().size(); j++) {

		SM += +poras[number].get_neighbor_Area(j);
		if (poras[poras[number].get_neighbor_Num(j)].get_filled()) {
			SMZ += poras[number].get_neighbor_Area(j);
		}
		if (!(poras[poras[number].get_neighbor_Num(j)].get_filled())) {
			surf_empty += poras[number].get_neighbor_Area(j);
		}
	}
	work = pres * 4.0 * pow(10.0, -9.0) * M_PI / 3.0 * pow(poras[number].get_Rad(), 3.0)
		- (4.0 * M_PI * poras[number].get_Rad() * poras[number].get_Rad() - SM) * DSIGMA
		- SIGMA * (surf_empty - SMZ);
	return work;
}
double rabotas_out_multi(std::vector<Pora>& poras, std::vector<int>& numbers, double pres, double DSIGMA, double SIGMA) {
	double SM = 0;
	double SMZ = 0;
	double surf_empty = 0;
	double nei_zap = 0;
	double work = 0;
	for (auto number : numbers) {
		for (int j = 0; j < poras[number].get_V_neighbor_Num().size(); j++) {

			SM += +poras[number].get_neighbor_Area(j);
			if (poras[poras[number].get_neighbor_Num(j)].get_filled()) {
				SMZ += poras[number].get_neighbor_Area(j);
			}
			if (!(poras[poras[number].get_neighbor_Num(j)].get_filled())) {
				surf_empty += poras[number].get_neighbor_Area(j);
			}
		}
		work += pres * 4.0 * pow(10.0, -9.0) * M_PI / 3.0 * pow(poras[number].get_Rad(), 3.0)
			- (4.0 * M_PI * poras[number].get_Rad() * poras[number].get_Rad() - SM) * DSIGMA
			- SIGMA * (surf_empty - SMZ);
	}

	for (auto number : numbers) {
		for (int j = 0; j < poras[number].get_V_neighbor_Num().size(); j++) {
			if (find_bool(numbers, poras[number].get_V_neighbor_Num()[j])) {
				work -= SIGMA * poras[number].get_neighbor_Area(j);
			}
		}
	}
	return work;

}

void filling(std::vector<Pora>& poras, std::vector<double>& volume_filled_graph, std::vector<double>& pressure_graph, double& volume_filled_sum, int N, int core_v, int DSIGMA_v, int persent_for_del) {
	double core = (1.0 * core_v) / 100.0;
	double DSIGMA = (1.0 * DSIGMA_v) / 1000.0;
	//std::string connectionString = "host=localhost port=1768 dbname=graphs user=postgres password =adu202121";
	int start_time = clock();
	double Rad = 5.0;     //средний радиус
	double sii = 1.0;   //дисперси€
	//double core = 0.8; // 0..1    1-касание, 0-все в одной точке
	//double DSIGMA = 40.0 / 1000.0; //DSIGMA поверхностна€ энерги€
	double SIGMA = 75.3 / 1000.0;//SIGMA поверхностна€ энерги€ жидкость-газ
	//int N = 30; //пор на ребре
	double L = 2.0 * Rad * core * N - Rad * core / 2.0;

	// параметры давлени€. Ќачало, конец, шаг на заполнение, шаг на вытекание
	double pressure_start = 0;
	double pressure_end = 34000001;
	double pressure_step = 100000;
	double poras_N = pow(N, 3.0);


	std::ifstream fin("fr_exp_9.dat");

	std::vector<double> x;
	std::vector<double> y;
	double d_tmp = 0;

	while (!fin.eof()) {

		fin >> d_tmp;
		x.push_back(d_tmp);
		fin >> d_tmp;
		y.push_back(d_tmp);
	}
	reverse(x.begin(), x.end());
	reverse(y.begin(), y.end());

	int N_theor = y.size();

	std::vector<double> interp_k;
	std::vector<double> interp_b;
	for (int i = 0; i < N_theor - 1; i++) {
		interp_k.push_back((y[i + 1] - y[i]) / (x[i + 1] - x[i]));
		interp_b.push_back(y[i] - interp_k[i] * x[i]);
	}



	/////////////////////////////////////////////////////////////////
	std::random_device rd_merge;
	std::mt19937 gen_merge(rd_merge());
	std::uniform_int_distribution<int> uid_merge(0, 99);
	int r_merge = 0;

	std::vector<Pora> poras_experiment;
	std::vector<Pora> poras_gauss;
	Generation_Poras_experiment(x, y, interp_k, interp_b, Rad, core, poras_experiment, N);
	Generation_Poras_gauss(Rad, core, poras_gauss, N, sii);
	
	for (int i = 0; i < poras_experiment.size(); i++) {
		r_merge = uid_merge(gen_merge);
		if (r_merge < 1000) {
			poras.push_back(poras_experiment[i]);
		}
		else {
			poras.push_back(poras_gauss[i]);
		}
	}
	/////////////////////////////////////////////////////////////////


	std::vector<int> bound_1;
	std::vector<int> bound_2;
	std::vector<int> bound_3;
	std::vector<int> bound_4;
	std::vector<int> bound_5;
	std::vector<int> bound_6;
	std::vector<int> bound;
	for (int i = 0; i < poras.size(); i++) {
		if (poras[i].get_X() == Rad * core) {
			poras[i].set_border(1);
			bound_1.push_back(i);
		}
		else if (poras[i].get_X() == Rad * core + (N - 1) * (Rad + Rad) * core) {
			poras[i].set_border(1);
			bound_2.push_back(i);
		}
		if (poras[i].get_Y() == Rad * core) {
			poras[i].set_border(1);
			bound_3.push_back(i);
		}
		else if (poras[i].get_Y() == Rad * core + (N - 1) * (Rad + Rad) * core) {
			poras[i].set_border(1);
			bound_4.push_back(i);
		}
		if (poras[i].get_Z() == Rad * core) {
			poras[i].set_border(1);
			bound_5.push_back(i);
		}
		else if (poras[i].get_Z() == Rad * core + (N - 1) * (Rad + Rad) * core) {
			poras[i].set_border(1);
			bound_6.push_back(i);
		}
		if (poras[i].get_border()) {
			bound.push_back(i);
		}

	}
	///////////////////////////////////////////////////////////////////

	
	std::vector<int> for_del(0);
	std::random_device rd_del;
	std::mt19937 gen_del(rd_del());
	std::uniform_int_distribution<int> uid_del(0, 99);
	int r_del = 0;
	for (int i = 0; i < poras.size(); i++) {
		if (!poras[i].get_border()) {
			r_del = uid_del(gen_del);
			if (r_del < persent_for_del) {
				for_del.push_back(i);
			}
		}
	}
	for (int i = 0; i < for_del.size(); i++) {
		poras[for_del[i]].set_Rad(0.001);
	}

	

	///////////////////////////////////////////////////////////////////


	int neighbors_pos[26] = { -1, +1, N + 1, N - 1, -N - 1, -N + 1, N, -N, N * N, N * N + 1, N * N - 1, N * N + N, N * N - N,
		   N * N + N + 1, N * N + N - 1, N * N - N + 1, N * N - N - 1, -N * N, -N * N + 1, -N * N - 1,
		   -N * N + N, -N * N - N, -N * N + N + 1, -N * N + N - 1, -N * N - N + 1, -N * N - N - 1 };
	int k;
	double kappa;
	double area_tmp;
	for (int i = 0; i < poras.size(); i += 2) {
		for (int j = 0; j < 26; j++) {
			k = i + neighbors_pos[j];
			if (poras[i].find_k(k)) {

			}
			else if (k <= (poras.size() - 1) && k >= 0) {
				kappa = pow((poras[i].get_X() - poras[k].get_X()) * (poras[i].get_X() - poras[k].get_X()) +
					(poras[i].get_Y() - poras[k].get_Y()) * (poras[i].get_Y() - poras[k].get_Y()) +
					(poras[i].get_Z() - poras[k].get_Z()) * (poras[i].get_Z() - poras[k].get_Z()), 0.5);
				if (kappa <= (poras[i].get_Rad() + poras[k].get_Rad())) {
					poras[i].set_neighbor_Num(k);
					poras[k].set_neighbor_Num(i);
					area_tmp = M_PI * ((-pow(poras[i].get_Rad(), 4.0) - pow(poras[k].get_Rad(), 4.0) - pow(kappa, 4.0) +
						2.0 * pow(poras[i].get_Rad(), 2.0) * pow(poras[k].get_Rad(), 2.0) +
						2.0 * pow(poras[i].get_Rad(), 2.0) * pow(kappa, 2.0) +
						2.0 * pow(kappa, 2.0) * pow(poras[k].get_Rad(), 2.0)) / (4.0 * kappa * kappa));
					poras[i].set_neighbor_Area(area_tmp);
					poras[k].set_neighbor_Area(area_tmp);
				}
				/*
				else
				{
					poras[i].set_neighbor_Num(-1);
					poras[i].set_neighbor_Area(0);
				}
				*/
			}
		}

		if ((i + 2) % (int)(N * N) == 0) {

		}
		else if ((i + 1) % (int)(N * N) == 0) {

		}
		else
			if ((i + 2) % (int)N == 0) {
				i = i + 1;
			}
			else if ((i + 1) % (int)N == 0) {
				i = i - 1;
			}
	}

	std::vector<int> zap_num;

	for (int i = 0; i < poras.size(); i++) {
		if (poras[i].get_border()) {
			poras[i].set_filled(1);
			zap_num.push_back(i);
		}
	}
	int zap_num_new = 0;
	double volume_filled = 0;
	int time_zap_start = clock();
	for (double pres = pressure_start; pres < pressure_end; pres += pressure_step) {
		zap_num_new = 1;
		//std::cout << pres << "\n";
		while (zap_num_new) {
			zap_num_new = 0;

			for (int j = 0; j < zap_num.size(); j++) {
				for (int k = 0; k < poras[zap_num[j]].get_V_neighbor_Num().size(); k++) {
					if (!(poras[poras[zap_num[j]].get_neighbor_Num(k)].get_filled())) {
						if (rabotas_in(poras, poras[zap_num[j]].get_neighbor_Num(k), pres, DSIGMA, SIGMA)) {
							poras[poras[zap_num[j]].get_neighbor_Num(k)].set_filled(1);
							volume_filled += ((4.0 * M_PI * pow(poras[poras[zap_num[j]].get_neighbor_Num(k)].get_Rad(), 3.0)) / 3.0);
							zap_num.push_back(poras[zap_num[j]].get_neighbor_Num(k));
							zap_num_new++;

						}

					}
				}
			}


		}
		volume_filled_graph.push_back(volume_filled);
		pressure_graph.push_back(pres / 101500.0);
	}

	int time_zap_end = clock();

	int time_zap_serch = time_zap_end - time_zap_start;
	//std::cout << N << "\n";
	//std::cout << time_zap_serch << "\n";
	double max_volume_filled;
	max_volume_filled = *max_element(volume_filled_graph.begin(), volume_filled_graph.end());
	volume_filled_sum = max_volume_filled;

	for (int i = 0; i < volume_filled_graph.size(); i++) {
		volume_filled_graph[i] = volume_filled_graph[i] / max_volume_filled;
	}


	int end_time = clock();
	int search_time = end_time - start_time;
	//std::cout << search_time << "\n";
}

void leakage(std::vector<Pora>& poras, std::vector<double>& volume_filled_graph_out, std::vector<double>& pressure_graph_out, double volume_filled, int core_v, int DSIGMA_v, int border_C, int outflow_M, std::vector<int>& porous_neibors) {
	double core = (1.0 * core_v) / 100.0;
	double DSIGMA = (1.0 * DSIGMA_v) / 1000.0;
	//std::string connectionString = "host=localhost port=1768 dbname=graphs user=postgres password =adu202121";
	int start_time = clock();
	std::vector<int> zap_num;
	for (int i = 0; i < poras.size(); i++) {
		if (poras[i].get_filled() && !(poras[i].get_border())) {
			zap_num.push_back(i);
		}
	}
	//////////////////////////////////////////////////////
	/*
	std::vector<int> for_del(0);
	std::random_device rd_del;
	std::mt19937 gen_del(rd_del());
	std::uniform_int_distribution<int> uid_del(0, 99);
	int r_del = 0;
	for (int i = 0; i < zap_num.size(); i++) {
		r_del = uid_del(gen_del);
		if (r_del < 3) {
			for_del.push_back(zap_num[i]);
		}
	}
	for (int i = 0; i < for_del.size(); i++) {
		volume_filled -= ((4.0 * M_PI * pow(poras[for_del[i]].get_Rad(), 3.0)) / 3.0);
		poras[for_del[i]].set_filled(0);
		auto it_del = std::find(zap_num.begin(), zap_num.end(), for_del[i]);
		int pos_zap_del = std::distance(begin(zap_num), it_del);
		zap_num.erase(zap_num.begin() + pos_zap_del);
	}
	*/
	//////////////////////////////////////////////////////
	double pressure_start = 0;
	double pressure_end = 34000001;
	double pressure_step = 100000;
	//double DSIGMA = 40.0 / 1000.0; //DSIGMA поверхностна€ энерги€
	double SIGMA = 75.3 / 1000.0;//SIGMA поверхностна€ энерги€ жидкость-газ
	double out_num_new = 0;
	std::vector<double> poras_Energy(0);
	std::vector<int> poras_Energy_pos(0);
	std::vector<double> poras_Energy_multi(0);
	std::vector<std::vector<int>> poras_Energy_pos_multi(0);
	double energy_out = 0.0;
	int pos_e = 0;
	int pos_e_multi = 0;
	int pos_e_p = 0;
	double e_level = 0.0;
	double e_level_multi = 0.0;
	for (double pres = pressure_end - 1; pres > pressure_start - 1; pres -= pressure_step) {
		/*
		if (outflow_M == 2) {
			std::cout << pres << std::endl;
		}
		*/
		out_num_new = 1;

		if (border_C == 1) {

			for (int j = 0; j < zap_num.size(); j++) {
				if (rabotas_out(poras, zap_num[j], pres, DSIGMA, SIGMA) <= 0.0) {
					if (poras[zap_num[j]].road_to_board(poras)) {
						for (int i = 0; i < poras.size(); i++) {
							poras[i].set_path_to_border(0);
						}
						poras[zap_num[j]].set_way_to_board(1);
					}
					else {
						for (int i = 0; i < poras.size(); i++) {
							poras[i].set_path_to_border(0);
						}
					}
				}
			}
			for (int j = 0; j < zap_num.size(); j++) {
				for (auto poraId : poras[zap_num[j]].get_V_neighbor_Num()) {
					if (poras[poraId].get_filled() && !poras[poraId].get_border()) {
						std::vector<int> tmp_vec = { zap_num[j], poraId };
						if (rabotas_out_multi(poras, tmp_vec, pres, DSIGMA, SIGMA) <= 0.0) {
							if (poras[zap_num[j]].road_to_board(poras)) {
								for (int i = 0; i < poras.size(); i++) {
									poras[i].set_path_to_border(0);
								}
								poras[zap_num[j]].set_way_to_board(1);
							}
							else {
								for (int i = 0; i < poras.size(); i++) {
									poras[i].set_path_to_border(0);
								}
							}

						}
					}
				}
			}
		}
		while (out_num_new) {
			out_num_new = 0;
			//for (int i = 1; i <= 26; i++) {
			poras_Energy.clear();
			poras_Energy.shrink_to_fit();
			poras_Energy_pos.clear();
			poras_Energy_pos.shrink_to_fit();
			poras_Energy_multi.clear();
			poras_Energy_multi.shrink_to_fit();
			poras_Energy_pos_multi.clear();
			poras_Energy_pos_multi.shrink_to_fit();
			for (int j = 0; j < zap_num.size(); j++) {
				int flag_way_to_board_pressure = 1;
				if (border_C == 1) {
					flag_way_to_board_pressure = poras[zap_num[j]].get_way_to_board();
				}
				if (flag_way_to_board_pressure) {
					if ((energy_out = rabotas_out(poras, zap_num[j], pres, DSIGMA, SIGMA)) <= 0.0) {
						int flag_way_to_board_each = 1;
						if (border_C == 2) {
							flag_way_to_board_each = poras[zap_num[j]].road_to_board(poras);
							for (int i = 0; i < poras.size(); i++) {
								poras[i].set_path_to_border(0);
							}
						}
						if (flag_way_to_board_each) {
							poras_Energy.push_back(energy_out);
							poras_Energy_pos.push_back(zap_num[j]);
						}
					}
				}
			}

			for (int j = 0; j < zap_num.size(); j++) {
				for (auto poraId : poras[zap_num[j]].get_V_neighbor_Num()) {
					if (poras[poraId].get_filled() && !poras[poraId].get_border()) {
						std::vector<int> tmp_vec = { zap_num[j], poraId };
						std::sort(tmp_vec.begin(), tmp_vec.end());
						if (!find_bool(poras_Energy_pos_multi, tmp_vec)) {
							int flag_way_to_board_pressure = 1;
							if (border_C == 1) {
								flag_way_to_board_pressure = poras[zap_num[j]].get_way_to_board();
							}
							if (flag_way_to_board_pressure) {
								if ((energy_out = rabotas_out_multi(poras, tmp_vec, pres, DSIGMA, SIGMA)) <= 0.0) {
									int flag_way_to_board_each = 1;
									if (border_C == 2) {
										flag_way_to_board_each = poras[zap_num[j]].road_to_board(poras);
										for (int i = 0; i < poras.size(); i++) {
											poras[i].set_path_to_border(0);
										}
									}
									if (flag_way_to_board_each) {
										poras_Energy_multi.push_back(energy_out);
										poras_Energy_pos_multi.push_back(tmp_vec);
									}
								}
							}
						}
					}
				}
			}

			poras_Energy_multi = {};
			if (poras_Energy.size() || poras_Energy_multi.size()) {
				switch (outflow_M) {
				case 0:
				{
					for (int k = 0; k < poras_Energy_pos.size(); k++) {
						pos_e_p = poras_Energy_pos[k];
						poras[pos_e_p].set_emptied(1);
						poras[pos_e_p].set_way_to_board(0);
						poras[pos_e_p].set_filled(0);
						volume_filled -= ((4.0 * M_PI * pow(poras[pos_e_p].get_Rad(), 3.0)) / 3.0);
						auto it = std::find(zap_num.begin(), zap_num.end(), pos_e_p);
						int pos_zap_num = std::distance(begin(zap_num), it);
						zap_num.erase(zap_num.begin() + pos_zap_num);
						out_num_new++;
						porous_neibors.push_back(poras[pos_e_p].filled_now(poras));
					}

					for (auto poras_Energy_pos : poras_Energy_pos_multi) {
						for (int k = 0; k < poras_Energy_pos.size(); k++) {
							pos_e_p = poras_Energy_pos[k];
							poras[pos_e_p].set_emptied(1);
							poras[pos_e_p].set_way_to_board(0);
							poras[pos_e_p].set_filled(0);
							auto it = std::find(zap_num.begin(), zap_num.end(), pos_e_p);
							if (it != std::end(zap_num)) {
								volume_filled -= ((4.0 * M_PI * pow(poras[pos_e_p].get_Rad(), 3.0)) / 3.0);
								int pos_zap_num = std::distance(begin(zap_num), it);
								zap_num.erase(zap_num.begin() + pos_zap_num);
								out_num_new++;
								porous_neibors.push_back(poras[pos_e_p].filled_now(poras));
							}
						}

					}
					break;
				}
				case 1:
				{
					std::random_device rd;
					std::mt19937 mt(rd());
					std::uniform_int_distribution<int> dist(0, poras_Energy.size() + poras_Energy_multi.size() - 1);
					pos_e = dist(mt);
					if (pos_e < poras_Energy.size()) {
						pos_e_p = poras_Energy_pos[pos_e];
						poras[pos_e_p].set_emptied(1);
						poras[pos_e_p].set_way_to_board(0);
						poras[pos_e_p].set_filled(0);
						volume_filled -= ((4.0 * M_PI * pow(poras[pos_e_p].get_Rad(), 3.0)) / 3.0);
						auto it = std::find(zap_num.begin(), zap_num.end(), pos_e_p);
						int pos_zap_num = std::distance(begin(zap_num), it);
						zap_num.erase(zap_num.begin() + pos_zap_num);
						out_num_new++;
						porous_neibors.push_back(poras[pos_e_p].filled_now(poras));
					}
					else {
						pos_e -= poras_Energy.size();
						for (auto pos_e_p : poras_Energy_pos_multi[pos_e]) {
							poras[pos_e_p].set_emptied(1);
							poras[pos_e_p].set_way_to_board(0);
							poras[pos_e_p].set_filled(0);
							volume_filled -= ((4.0 * M_PI * pow(poras[pos_e_p].get_Rad(), 3.0)) / 3.0);
							auto it = std::find(zap_num.begin(), zap_num.end(), pos_e_p);
							int pos_zap_num = std::distance(begin(zap_num), it);
							zap_num.erase(zap_num.begin() + pos_zap_num);
							out_num_new++;
							porous_neibors.push_back(poras[pos_e_p].filled_now(poras));
						}
					}
					//std::vector<double>::iterator result = select_randomly(poras_Energy.begin(), poras_Energy.end());
					//pos_e = std::distance(begin(poras_Energy), result);
					

					break;
				}
				case 2:
				{
					std::vector<double>::iterator result;
					std::vector<double>::iterator result_multi;
					if (poras_Energy.size()) {
						result = std::min_element(begin(poras_Energy), end(poras_Energy));
						pos_e = std::distance(begin(poras_Energy), result);
						e_level = poras_Energy[pos_e];
					}
					else {
						e_level = 100;
					}
					if (poras_Energy_multi.size()) {
						result_multi = std::min_element(begin(poras_Energy_multi), end(poras_Energy_multi));
						pos_e_multi = std::distance(begin(poras_Energy_multi), result_multi);
						e_level_multi = poras_Energy_multi[pos_e_multi];
					}
					else {
						e_level_multi = 100;
					}

					if (e_level < e_level_multi) {
						pos_e_p = poras_Energy_pos[pos_e];
						poras[pos_e_p].set_emptied(1);
						poras[pos_e_p].set_way_to_board(0);
						poras[pos_e_p].set_filled(0);
						volume_filled -= ((4.0 * M_PI * pow(poras[pos_e_p].get_Rad(), 3.0)) / 3.0);
						auto it = std::find(zap_num.begin(), zap_num.end(), pos_e_p);
						int pos_zap_num = std::distance(zap_num.begin(), it);
						zap_num.erase(zap_num.begin() + pos_zap_num);
						out_num_new++;
						porous_neibors.push_back(poras[pos_e_p].filled_now(poras));
					}
					else {
						//std::cout << "e_level_multi: "<< e_level_multi<<" e_level: "<< e_level << " pres: " << pres << " zap_num: "<<zap_num.size() << std::endl;
						for (auto pos_e_p : poras_Energy_pos_multi[pos_e_multi]) {
							//summmm++;
							poras[pos_e_p].set_emptied(1);
							poras[pos_e_p].set_way_to_board(0);
							poras[pos_e_p].set_filled(0);
							volume_filled -= ((4.0 * M_PI * pow(poras[pos_e_p].get_Rad(), 3.0)) / 3.0);
							auto it = std::find(zap_num.begin(), zap_num.end(), pos_e_p);
							int pos_zap_num = std::distance(zap_num.begin(), it);
							/*
							if (pres <= 4100002) {
								std::cout << pos_e_p << std::endl;
								std::cout << pos_zap_num << std::endl;
								std::cout << "check" << std::endl;
								
								std::ofstream fout;
								fout.open("txt.txt");
								for (auto num : zap_num) {
									fout << num<< std::endl;
								}
								fout.close();
								
							}
							*/
							zap_num.erase(zap_num.begin() + pos_zap_num);
							out_num_new++;
							porous_neibors.push_back(poras[pos_e_p].filled_now(poras));
						}
					}
					break;
				}
				default:
					break;
				}
				if (border_C == 2) {
					for (int i = 0; i < zap_num.size(); i++) {
						poras[zap_num[i]].set_way_to_board(0);
					}
				}
				/*
				//std::vector<double>::iterator result = select_randomly(poras_Energy.begin(), poras_Energy.end());
				std::vector<double>::iterator result = std::min_element(begin(poras_Energy), end(poras_Energy));

				pos_e = std::distance(begin(poras_Energy), result);
				pos_e_p = poras_Energy_pos[pos_e];
				poras[pos_e_p].set_emptied(1);
				poras[pos_e_p].set_filled(0);
				volume_filled -= ((4.0 * M_PI * pow(poras[pos_e_p].get_Rad(), 3.0)) / 3.0);
				auto it = std::find(zap_num.begin(), zap_num.end(), pos_e_p);
				int pos_zap_num = std::distance(begin(zap_num), it);
				zap_num.erase(zap_num.begin() + pos_zap_num);
				//j = -1;
				//i = 0;
				out_num_new++;
				*/

			}
			//}

		}
		for (int i = 0; i < zap_num.size(); i++) {
			poras[zap_num[i]].set_way_to_board(0);
		}
		volume_filled_graph_out.push_back(volume_filled);
		pressure_graph_out.push_back(pres / 101500.0);
	}
	double max_volume_filled;
	max_volume_filled = *max_element(volume_filled_graph_out.begin(), volume_filled_graph_out.end());
	for (int i = 0; i < volume_filled_graph_out.size(); i++) {
		volume_filled_graph_out[i] = volume_filled_graph_out[i] / max_volume_filled;
	}
	int end_time = clock();
	int search_time = end_time - start_time;
	//std::cout << search_time << "\n";
	/*
	std::cout << porous_neibors.size() << std::endl;
	for (auto object : porous_neibors)
	{
		std::cout << object << " ";
	}
	std::cout << std::endl;
	*/
}

void genegate_body(int N, std::vector<int> core_v, std::vector<int> DSIGMA_v, int i_min, int i_max, std::string connectionString, int body_id, std::vector<int> persent_for_del_vec) {
	mtx.lock();
	std::cout << "Generate id: " << std::this_thread::get_id() << " started" << std::endl;
	mtx.unlock();
	for (size_t i = i_min; i < i_max; i++) {
		for (size_t j = 0; j < DSIGMA_v.size(); j++) {
			for (auto persent_for_del : persent_for_del_vec) {
				std::string table_name = "body_core_" + std::to_string(core_v[i]) + "_DS_" + std::to_string(DSIGMA_v[j]) + "_ID_" + std::to_string(body_id)+ "_PFD_"+ std::to_string(persent_for_del);
				std::string graph_id = std::to_string(core_v[i]) + std::to_string(DSIGMA_v[j]) + std::to_string(body_id) + std::to_string(persent_for_del);
				std::vector<Pora> poras;
				std::vector<double> volume_filled_graph;
				std::vector<double> pressure_graph;
				double volume_filled = 0;
				filling(poras, volume_filled_graph, pressure_graph, volume_filled, N, core_v[i], DSIGMA_v[j], persent_for_del);
				std::string connectionString_graphs = "host=localhost port=1768 dbname=graphs user=postgres password =adu202121";
				insert_into_table_graph(graph_id, volume_filled_graph, pressure_graph, "graph_in", connectionString_graphs);
				create_table(table_name, connectionString);
				insert_into_table(poras, table_name, connectionString);

			}
		}
	}
	mtx.lock();
	std::cout << "Generate id: " << std::this_thread::get_id() << " finished" << std::endl;
	mtx.unlock();
}

void empty_body(std::vector<int> core_v, std::vector<int> DSIGMA_v, int i_min, int i_max, std::string connectionString, std::vector<int> body_id, int border_C, int outflow_M, std::vector<int> persent_for_del_vec) {
	mtx.lock();
	std::cout << "Empty id: " << std::this_thread::get_id() << " " << border_C << " " << outflow_M << " started" << std::endl;
	mtx.unlock();
	for (size_t k = 0; k < body_id.size(); k++) {
		for (size_t i = i_min; i < i_max; i++) {
			for (size_t j = 0; j < DSIGMA_v.size(); j++) {
				for (auto persent_for_del : persent_for_del_vec) {
					std::string table_name = "body_core_" + std::to_string(core_v[i]) + "_DS_" + std::to_string(DSIGMA_v[j]) + "_ID_" + std::to_string(body_id[k]) + "_PFD_" + std::to_string(persent_for_del);
					std::string graph_id = std::to_string(core_v[i]) + std::to_string(DSIGMA_v[j]) + std::to_string(body_id[k]) + std::to_string(border_C) + std::to_string(outflow_M) + std::to_string(persent_for_del);
					std::vector<Pora> poras;
					select_from_table(poras, table_name, connectionString);
					double volume_filled = 0;

					for (size_t i = 0; i < poras.size(); i++) {
						if (poras[i].get_filled() && !(poras[i].get_border())) {
							volume_filled += ((4.0 * M_PI * pow(poras[i].get_Rad(), 3.0)) / 3.0);
						}
					}
					std::vector<double> volume_filled_graph_out;
					std::vector<double> pressure_graph_out;
					std::vector<int> porous_neibors;
					leakage(poras, volume_filled_graph_out, pressure_graph_out, volume_filled, core_v[i], DSIGMA_v[j], border_C, outflow_M, porous_neibors);

					std::string connectionString_graphs = "host=localhost port=1768 dbname=graphs user=postgres password =adu202121";
					insert_into_table_graph(graph_id, volume_filled_graph_out, pressure_graph_out, "graph_out", connectionString_graphs);
					connectionString_graphs = "host=localhost port=1768 dbname=graphs user=postgres password =adu202121";
					insert_into_table_graph_li(graph_id, porous_neibors, connectionString_graphs);
				}
			}
		}
	}
	mtx.lock();
	std::cout << "Empty id: " << std::this_thread::get_id() << " " << border_C << " " << outflow_M << " finished" << std::endl;;
	mtx.unlock();
}


int main() {
	///
	///
	/// ѕовысить core что позволит уменьшить количство сеседей, но сильно уменьшить DSigma
	/// 
	/// ѕосле генерации тела по гаусовому закону, пройти по порам и у пор с большим количеством соседей заменить радиус случайной величининой из  эксперимента
	/// 
	///
	///
	/// 
	/// 
	///
	SimpleTimer st;
	SetConsoleOutputCP(1251);
	std::string connectionString = "host=localhost port=1768 dbname=porous_body user=postgres password =adu202121";

	//pqxx::result response;
	/*
	pqxx::result response = worker.exec("SELECT * FROM body_1 WHERE id=3");

	for (size_t i = 0; i < response.size(); i++)
	{
		for (size_t j = 0; j < response[i].size(); j++)
		{
			str += to_string(response[i][j])+" ";
		}

		std::cout <<str<< std::endl;

	}
*/

// border_C: 0 - не учитываем, 1 - расчЄт после измениени€ давлени€, 2 - расчЄт дл€ каждой
// outflow_M: 0 - все, 1 - случайна€, 2 - с минимальной энерги€ей
//std::vector<int> body_id = { 1,2,3,4,5,6,7,8,9,10 };
	std::vector<int> body_id = { 66 };
	std::vector<int> border_C = { 0,1,2 };
	std::vector<int> outflow_M = { 0,1,2 };
	std::vector<int> core_v;
	std::vector<int> DSIGMA_v;
	for (size_t i = 0; i <= 10; i++) {
		core_v.push_back(93 + 2 * i);
	}
	for (size_t i = 0; i <= 0; i++) {
		DSIGMA_v.push_back(9 + 2 * i);
	}

	std::vector<int> persent_for_del_vec = { 0 };

	std::cout << "ѕор на границе: ";
	int N;
	std::cin >> N;

	/////////////////////////////////////////////////
	std::thread th0_in(genegate_body, N, core_v, DSIGMA_v, 0 * core_v.size() / 11, 1 * core_v.size() / 11, connectionString, body_id[0], persent_for_del_vec);
	
	//std::thread th1_in(genegate_body, N, core_v, DSIGMA_v, 1 * core_v.size() / 11, 2 * core_v.size() / 11, connectionString, body_id[0], persent_for_del_vec);
	/*
	std::thread th2_in(genegate_body, N, core_v, DSIGMA_v, 2 * core_v.size() / 11, 3 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th3_in(genegate_body, N, core_v, DSIGMA_v, 3 * core_v.size() / 11, 4 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th4_in(genegate_body, N, core_v, DSIGMA_v, 4 * core_v.size() / 11, 5 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th5_in(genegate_body, N, core_v, DSIGMA_v, 5 * core_v.size() / 11, 6 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th6_in(genegate_body, N, core_v, DSIGMA_v, 6 * core_v.size() / 11, 7 * core_v.size() / 11, connectionString, body_id[0]);
	*/
	/*
	std::thread th7_in(genegate_body, N, core_v, DSIGMA_v, 7 * core_v.size() / 11, 8 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th8_in(genegate_body, N, core_v, DSIGMA_v, 8 * core_v.size() / 11, 9 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th9_in(genegate_body, N, core_v, DSIGMA_v, 9 * core_v.size() / 11, 10 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th10_in(genegate_body, N, core_v, DSIGMA_v, 10 * core_v.size() / 11, 11 * core_v.size() / 11, connectionString, body_id[0]);
	*/

	/////////////////////////////////////////////////
	th0_in.join();
	
	//th1_in.join();
	/*
	th2_in.join();
	th3_in.join();
	th4_in.join();
	th5_in.join();
	th6_in.join();
	*/
	/*
	th7_in.join();
	th8_in.join();
	th9_in.join();
	th10_in.join();
	*/

	/////////////////////////////////////////////////

	//////////////////////////////////////////
	
	std::thread th00_out(empty_body, core_v, DSIGMA_v, 0 * core_v.size() / 11, 1 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0], persent_for_del_vec);
	/*
	std::thread th01_out(empty_body, core_v, DSIGMA_v, 1 * core_v.size() / 11, 2 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th02_out(empty_body, core_v, DSIGMA_v, 2 * core_v.size() / 11, 3 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th03_out(empty_body, core_v, DSIGMA_v, 3 * core_v.size() / 11, 4 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th04_out(empty_body, core_v, DSIGMA_v, 4 * core_v.size() / 11, 5 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th05_out(empty_body, core_v, DSIGMA_v, 5 * core_v.size() / 11, 6 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th06_out(empty_body, core_v, DSIGMA_v, 6 * core_v.size() / 11, 7 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	*/
	/*
	std::thread th07_out(empty_body, core_v, DSIGMA_v, 7 * core_v.size() / 11, 8 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th08_out(empty_body, core_v, DSIGMA_v, 8 * core_v.size() / 11, 9 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th09_out(empty_body, core_v, DSIGMA_v, 9 * core_v.size() / 11, 10 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th010_out(empty_body, core_v, DSIGMA_v, 10 * core_v.size() / 11, 11 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	*/
	
	//////////////////////////////////////////
	
	//th00_out.join();
	/*
	th01_out.join();
	th02_out.join();
	th03_out.join();
	th04_out.join();
	th05_out.join();
	th06_out.join();
	*/
	/*
	th07_out.join();
	th08_out.join();
	th09_out.join();
	th010_out.join();
	*/
	
	//////////////////////////////////////////
	
	std::thread th10_out(empty_body, core_v, DSIGMA_v, 0 * core_v.size() / 11, 1 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1], persent_for_del_vec);
	/*
	std::thread th11_out(empty_body, core_v, DSIGMA_v, 1 * core_v.size() / 11, 2 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th12_out(empty_body, core_v, DSIGMA_v, 2 * core_v.size() / 11, 3 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th13_out(empty_body, core_v, DSIGMA_v, 3 * core_v.size() / 11, 4 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th14_out(empty_body, core_v, DSIGMA_v, 4 * core_v.size() / 11, 5 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th15_out(empty_body, core_v, DSIGMA_v, 5 * core_v.size() / 11, 6 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th16_out(empty_body, core_v, DSIGMA_v, 6 * core_v.size() / 11, 7 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	*/
	/*
	std::thread th17_out(empty_body, core_v, DSIGMA_v, 7 * core_v.size() / 11, 8 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th18_out(empty_body, core_v, DSIGMA_v, 8 * core_v.size() / 11, 9 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th19_out(empty_body, core_v, DSIGMA_v, 9 * core_v.size() / 11, 10 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th110_out(empty_body, core_v, DSIGMA_v, 10 * core_v.size() / 11, 11 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	*/

	//////////////////////////////////////////
	
	//th10_out.join();
	/*
	th11_out.join();
	th12_out.join();
	th13_out.join();
	th14_out.join();
	th15_out.join();
	th16_out.join();
	*/
	/*
	th17_out.join();
	th18_out.join();
	th19_out.join();
	th110_out.join();
	*/

	//////////////////////////////////////////
	
	std::thread th20_out(empty_body, core_v, DSIGMA_v, 0 * core_v.size() / 11, 1 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2], persent_for_del_vec);
	/*
	std::thread th21_out(empty_body, core_v, DSIGMA_v, 1 * core_v.size() / 11, 2 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th22_out(empty_body, core_v, DSIGMA_v, 2 * core_v.size() / 11, 3 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th23_out(empty_body, core_v, DSIGMA_v, 3 * core_v.size() / 11, 4 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th24_out(empty_body, core_v, DSIGMA_v, 4 * core_v.size() / 11, 5 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th25_out(empty_body, core_v, DSIGMA_v, 5 * core_v.size() / 11, 6 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th26_out(empty_body, core_v, DSIGMA_v, 6 * core_v.size() / 11, 7 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	*/
	/*
	std::thread th27_out(empty_body, core_v, DSIGMA_v, 7 * core_v.size() / 11, 8 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th28_out(empty_body, core_v, DSIGMA_v, 8 * core_v.size() / 11, 9 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th29_out(empty_body, core_v, DSIGMA_v, 9 * core_v.size() / 11, 10 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th210_out(empty_body, core_v, DSIGMA_v, 10 * core_v.size() / 11, 11 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	*/

	//////////////////////////////////////////
	th00_out.join();
	th10_out.join();
	th20_out.join();
	/*
	th21_out.join();
	th22_out.join();
	th23_out.join();
	th24_out.join();
	th25_out.join();
	th26_out.join();
	*/
	/*
	th27_out.join();
	th28_out.join();
	th29_out.join();
	th210_out.join();
	*/

	//////////////////////////////////////////
	
	//////////////////////////////////////////
	//////////////////////////////////////////

	std::cout << "Success" << std::endl;

	return 0;
}
