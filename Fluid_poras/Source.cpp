
# define M_PI           3.14159265358979323846  /* pi */

#pragma comment(linker, "/STACK:128000000")

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
//#include <ctime>
#include <algorithm>
#include <iterator>
#include <pqxx/pqxx>
#include <Windows.h>
#include <thread>
#include <mutex>
#include "pore.h"
#include "supportFunc.h"
#include "pSql.h"
#include "SimpleTimer.h"

std::mutex mtx;
//using namespace std;
//int summmm = 0;

void generation_pores_gauss(const double rad, const double core, std::vector<pore>& pores, const double n, const double sii) {
	std::random_device rd_r;
	std::normal_distribution<> uid_r(rad, sii);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     uid_R(rd_R)
				pores.emplace_back(rad * core + 1.0 * i * (rad + rad) * core, rad * core + 1.0 * j * (rad + rad) * core, rad * core + 1.0 * k * (rad + rad) * core, uid_r(rd_r), 0, 0, 0, 0, 0);
			}
		}
	}
}
void generation_pores_experiment(std::vector<double>& x, std::vector<double>& y, std::vector<double>& interp_k, std::vector<double>& interp_b, double Rad, double core, std::vector<pore>& pores, double N) {
	const double min_x = *std::ranges::min_element(x);
	const double max_x = *std::ranges::max_element(x);
	const double min_y = *std::ranges::min_element(y);
	const double max_y = *std::ranges::max_element(y);
	std::random_device rd_x;
	std::mt19937 mt_x(rd_x());
	const std::uniform_real_distribution<> dist_x(min_x, max_x);
	//pos_e = dist(mt);
	std::random_device rd_y;
	std::mt19937 mt_y(rd_y());
	const std::uniform_real_distribution<> dist_y(min_y, max_y);
	double rad_tmp = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				bool kj = false;
				while (!kj) {
					if (const double pot_por_rad = dist_x(mt_x); dist_y(mt_y) <= linear(x, y, interp_k, interp_b, pot_por_rad, x.size())) {
						rad_tmp = pot_por_rad;
						kj = true;
					}
				}
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     Rad_tmp
				pores.emplace_back(Rad * core + 1.0 * i * (Rad + Rad) * core, Rad * core + 1.0 * j * (Rad + Rad) * core, Rad * core + 1.0 * k * (Rad + Rad) * core, rad_tmp, 0, 0, 0, 0, 0);
			}
		}
	}
}

int rabotas_in(std::vector<pore>& pores, const int i, const double pres, const double d_sigma, const double sigma) {
	double sm = 0;
	double smz = 0;
	double surf_empty = 0;
	double nei_zap = 0;
	
	for (size_t j = 0; j < pores[i].get_v_neighbor_num().size(); j++) {

		sm += +pores[i].get_neighbor_area(j);
		if (pores[pores[i].get_neighbor_num(j)].get_filled()) {
			nei_zap++;
			smz += pores[i].get_neighbor_area(j);
		}
		if (!(pores[pores[i].get_neighbor_num(j)].get_filled())) {
			surf_empty += pores[i].get_neighbor_area(j);
		}
	}
	const double work = -pres * 4.0 * pow(10.0, -9.0) * M_PI / 3.0 * pow(pores[i].get_rad(), 3.0) +
		(4.0 * M_PI * pores[i].get_rad() * pores[i].get_rad() - sm) * d_sigma +
		sigma * (surf_empty - smz);
	return (work < 0);
}

double rabotas_out(std::vector<pore>& pores, const int number, const double pres, const double d_sigma, const double sigma) {
	double sm = 0;
	double smz = 0;
	double surf_empty = 0;
	
	for (size_t j = 0; j < pores[number].get_v_neighbor_num().size(); j++) {

		sm += +pores[number].get_neighbor_area(j);
		if (pores[pores[number].get_neighbor_num(j)].get_filled()) {
			smz += pores[number].get_neighbor_area(j);
		}
		if (!(pores[pores[number].get_neighbor_num(j)].get_filled())) {
			surf_empty += pores[number].get_neighbor_area(j);
		}
	}
	const double work = pres * 4.0 * pow(10.0, -9.0) * M_PI / 3.0 * pow(pores[number].get_rad(), 3.0)
		- (4.0 * M_PI * pores[number].get_rad() * pores[number].get_rad() - sm) * d_sigma
		- sigma * (surf_empty - smz);
	return work;
}
double rabotas_out_multi(std::vector<pore>& pores, std::vector<int>& numbers, const double pres, const double d_sigma, const double sigma) {
	double sm = 0;
	double smz = 0;
	double surf_empty = 0;
	double work = 0;
	for (const auto number : numbers) {
		for (size_t j = 0; j < pores[number].get_v_neighbor_num().size(); j++) {

			sm += +pores[number].get_neighbor_area(j);
			if (pores[pores[number].get_neighbor_num(j)].get_filled()) {
				smz += pores[number].get_neighbor_area(j);
			}
			if (!(pores[pores[number].get_neighbor_num(j)].get_filled())) {
				surf_empty += pores[number].get_neighbor_area(j);
			}
		}
		work += pres * 4.0 * pow(10.0, -9.0) * M_PI / 3.0 * pow(pores[number].get_rad(), 3.0)
			- (4.0 * M_PI * pores[number].get_rad() * pores[number].get_rad() - sm) * d_sigma
			- sigma * (surf_empty - smz);
	}

	for (const auto number : numbers) {
		for (size_t j = 0; j < pores[number].get_v_neighbor_num().size(); j++) {
			if (find_bool(numbers, pores[number].get_v_neighbor_num()[j])) {
				work -= sigma * pores[number].get_neighbor_area(j);
			}
		}
	}
	return work;

}

void filling(std::vector<pore>& pores, std::vector<double>& volume_filled_graph, std::vector<double>& pressure_graph, double& volume_filled_sum, int n, int core_v, int d_sigma_v, int percent_for_del, std::vector<int>& porous_neighbours) {
	double core = (1.0 * core_v) / 100.0;
	double d_sigma = (1.0 * d_sigma_v) / 1000.0;
	
	double rad = 5.0;     //средний радиус
	double sii = 1.0;   //дисперси€
	//double core = 0.8; // 0..1    1-касание, 0-все в одной точке
	//double d_sigma = 40.0 / 1000.0; //d_sigma поверхностна€ энерги€
	double sigma = 74.23 / 1000.0;//SIGMA поверхностна€ энерги€ жидкость-газ
	//int N = 30; //пор на ребре
	//double L = 2.0 * rad * core * N - rad * core / 2.0;

	// параметры давлени€. Ќачало, конец, шаг на заполнение, шаг на вытекание
	double pressure_start = 0;
	double pressure_end = 34000001;
	double pressure_step = 100000;
	//double pores_N = pow(N, 3.0);


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
	std::ranges::reverse(x);
	std::ranges::reverse(y);

	size_t n_theor = y.size();

	std::vector<double> interp_k;
	std::vector<double> interp_b;
	for (size_t i = 0; i < n_theor - 1; i++) {
		interp_k.push_back((y[i + 1] - y[i]) / (x[i + 1] - x[i]));
		interp_b.push_back(y[i] - interp_k[i] * x[i]);
	}



	/////////////////////////////////////////////////////////////////
	std::random_device rd_merge;
	std::mt19937 gen_merge(rd_merge());
	std::uniform_int_distribution<int> uid_merge(0, 99);
	int r_merge = 0;

	std::vector<pore> pores_experiment;
	std::vector<pore> pores_gauss;
	generation_pores_experiment(x, y, interp_k, interp_b, rad, core, pores_experiment, n);
	generation_pores_gauss(rad, core, pores_gauss, n, sii);
	
	for (size_t i = 0; i < pores_experiment.size(); i++) {
		r_merge = uid_merge(gen_merge);
		if (r_merge < -1000) {
			pores.push_back(pores_experiment[i]);
		}
		else {
			pores.push_back(pores_gauss[i]);
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
	for (size_t i = 0; i < pores.size(); i++) {
		if (pores[i].get_x() == rad * core) {
			pores[i].set_border(1);
			bound_1.push_back(i);
		}
		else if (pores[i].get_x() == rad * core + (n - 1) * (rad + rad) * core) {
			pores[i].set_border(1);
			bound_2.push_back(i);
		}
		if (pores[i].get_y() == rad * core) {
			pores[i].set_border(1);
			bound_3.push_back(i);
		}
		else if (pores[i].get_y() == rad * core + (n - 1) * (rad + rad) * core) {
			pores[i].set_border(1);
			bound_4.push_back(i);
		}
		if (pores[i].get_z() == rad * core) {
			pores[i].set_border(1);
			bound_5.push_back(i);
		}
		else if (pores[i].get_z() == rad * core + (n - 1) * (rad + rad) * core) {
			pores[i].set_border(1);
			bound_6.push_back(i);
		}
		if (pores[i].get_border()) {
			bound.push_back(i);
		}

	}
	///////////////////////////////////////////////////////////////////

	
	std::vector<int> for_del(0);
	std::random_device rd_del;
	std::mt19937 gen_del(rd_del());
	std::uniform_int_distribution<int> uid_del(0, 99);
	int r_del = 0;
	for (size_t i = 0; i < pores.size(); i++) {
		if (true  /* && !pores[i].get_border()*/) {
			r_del = uid_del(gen_del);
			if (r_del < percent_for_del) {
				for_del.push_back(i);
			}
		}
	}
	for (int i : for_del)
	{
		pores[i].set_rad(0.001);
	}

	

	///////////////////////////////////////////////////////////////////


	int neighbors_pos[26] = { -1, +1, n + 1, n - 1, -n - 1, -n + 1, n, -n, n * n, n * n + 1, n * n - 1, n * n + n, n * n - n,
		   n * n + n + 1, n * n + n - 1, n * n - n + 1, n * n - n - 1, -n * n, -n * n + 1, -n * n - 1,
		   -n * n + n, -n * n - n, -n * n + n + 1, -n * n + n - 1, -n * n - n + 1, -n * n - n - 1 };
	int k;
	double kappa;
	double area_tmp;
	for (size_t i = 0; i < pores.size(); i += 2) {
		for (int neighbors_po : neighbors_pos)
		{
			k = i + neighbors_po;
			if (pores[i].find_k(k)) {

			}
			else if (k <= (pores.size() - 1) && k >= 0) {
				kappa = pow((pores[i].get_x() - pores[k].get_x()) * (pores[i].get_x() - pores[k].get_x()) +
					(pores[i].get_y() - pores[k].get_y()) * (pores[i].get_y() - pores[k].get_y()) +
					(pores[i].get_z() - pores[k].get_z()) * (pores[i].get_z() - pores[k].get_z()), 0.5);
				if (kappa <= (pores[i].get_rad() + pores[k].get_rad())) {
					pores[i].set_neighbor_num(k);
					pores[k].set_neighbor_num(i);
					area_tmp = M_PI * ((-pow(pores[i].get_rad(), 4.0) - pow(pores[k].get_rad(), 4.0) - pow(kappa, 4.0) +
						2.0 * pow(pores[i].get_rad(), 2.0) * pow(pores[k].get_rad(), 2.0) +
						2.0 * pow(pores[i].get_rad(), 2.0) * pow(kappa, 2.0) +
						2.0 * pow(kappa, 2.0) * pow(pores[k].get_rad(), 2.0)) / (4.0 * kappa * kappa));
					pores[i].set_neighbor_area(area_tmp);
					pores[k].set_neighbor_area(area_tmp);
				}
				/*
				else
				{
					pores[i].set_neighbor_Num(-1);
					pores[i].set_neighbor_Area(0);
				}
				*/
			}
		}

		if ((i + 2) % (static_cast<unsigned long long>(n) * n) == 0) {

		}
		else if ((i + 1) % (static_cast<unsigned long long>(n) * n) == 0) {

		}
		else
			if ((i + 2) % n == 0) {
				i = i + 1;
			}
			else if ((i + 1) % n == 0) {
				i = i - 1;
			}
	}

	std::vector<int> zap_num;

	for (size_t i = 0; i < pores.size(); i++) {
		if (pores[i].get_border()) {
			pores[i].set_filled(1);
			zap_num.push_back(i);
		}
	}
	int zap_num_new = 0;
	double volume_filled = 0;
	
	for (double pres = pressure_start; pres < pressure_end; pres += pressure_step) {
		zap_num_new = 1;
		//std::cout << pres << "\n";
		while (zap_num_new) {
			zap_num_new = 0;

			for (size_t j = 0; j < zap_num.size(); j++) {
				for (size_t size = 0; size < pores[zap_num[j]].get_v_neighbor_num().size(); size++) {
					if (!(pores[pores[zap_num[j]].get_neighbor_num(size)].get_filled())) {
						if (rabotas_in(pores, pores[zap_num[j]].get_neighbor_num(size), pres, d_sigma, sigma)) {
							pores[pores[zap_num[j]].get_neighbor_num(size)].set_filled(1);
							volume_filled += ((4.0 * M_PI * pow(pores[pores[zap_num[j]].get_neighbor_num(size)].get_rad(), 3.0)) / 3.0);
							zap_num.push_back(pores[zap_num[j]].get_neighbor_num(size));
							zap_num_new++;
							porous_neighbours.push_back(pores[pores[zap_num[j]].get_neighbor_num(size)].filled_now(pores));

						}

					}
				}
			}


		}
		volume_filled_graph.push_back(volume_filled);
		pressure_graph.push_back(pres / 101500.0);
	}

	


	//std::cout << N << "\n";
	//std::cout << time_zap_serch << "\n";
	double max_volume_filled;
	max_volume_filled = *std::ranges::max_element(volume_filled_graph);
	volume_filled_sum = max_volume_filled;

	for (double& i : volume_filled_graph)
	{
		i = i / max_volume_filled;
	}


	

	//std::cout << search_time << "\n";
}

void leakage(std::vector<pore>& pores, std::vector<double>& volume_filled_graph_out,
             std::vector<double>& pressure_graph_out, double volume_filled, int core_v, int d_sigma_v, int border_c,
             int outflow_m, std::vector<int>& porous_neighbours)
{
	
	double d_sigma = (1.0 * d_sigma_v) / 1000.0;
	
	std::vector<int> zap_num;
	for (size_t i = 0; i < pores.size(); i++) {
		if (pores[i].get_filled() && !(pores[i].get_border())) {
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
		volume_filled -= ((4.0 * M_PI * pow(pores[for_del[i]].get_Rad(), 3.0)) / 3.0);
		pores[for_del[i]].set_filled(0);
		auto it_del = std::find(zap_num.begin(), zap_num.end(), for_del[i]);
		int pos_zap_del = std::distance(begin(zap_num), it_del);
		zap_num.erase(zap_num.begin() + pos_zap_del);
	}
	*/
	//////////////////////////////////////////////////////
	double pressure_start = 0;
	double pressure_end = 34000001;
	double pressure_step = 100000;
	//double d_sigma = 40.0 / 1000.0; //d_sigma поверхностна€ энерги€
	double sigma = 74.23 / 1000.0;//SIGMA поверхностна€ энерги€ жидкость-газ
	int out_num_new = 0;
	std::vector<double> pores_energy(0);
	std::vector<int> pores_energy_pos(0);
	std::vector<double> pores_energy_multi(0);
	std::vector<std::vector<int>> pores_energy_pos_multi(0);
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

		if (border_c == 1) {

			for (int j : zap_num)
			{
				if (rabotas_out(pores, j, pres, d_sigma, sigma) <= 0.0) {
					if (pores[j].road_to_board(pores)) {
						for (auto& pore : pores)
						{
							pore.set_path_to_border(0);
						}
						pores[j].set_way_to_board(1);
					}
					else {
						for (auto& pore : pores)
						{
							pore.set_path_to_border(0);
						}
					}
				}
			}
			for (int& j : zap_num)
			{
				for (auto poreId : pores[j].get_v_neighbor_num()) {
					if (pores[poreId].get_filled() && !pores[poreId].get_border()) {
						std::vector<int> tmp_vec = {j, poreId };
						if (rabotas_out_multi(pores, tmp_vec, pres, d_sigma, sigma) <= 0.0) {
							if (pores[j].road_to_board(pores)) {
								for (auto& pore : pores)
								{
									pore.set_path_to_border(0);
								}
								pores[j].set_way_to_board(1);
							}
							else {
								for (auto& pore : pores)
								{
									pore.set_path_to_border(0);
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
			pores_energy.clear();
			pores_energy.shrink_to_fit();
			pores_energy_pos.clear();
			pores_energy_pos.shrink_to_fit();
			pores_energy_multi.clear();
			pores_energy_multi.shrink_to_fit();
			pores_energy_pos_multi.clear();
			pores_energy_pos_multi.shrink_to_fit();
			for (int j : zap_num)
			{
				int flag_way_to_board_pressure = 1;
				if (border_c == 1) {
					flag_way_to_board_pressure = pores[j].get_way_to_board();
				}
				if (flag_way_to_board_pressure) {
					if ((energy_out = rabotas_out(pores, j, pres, d_sigma, sigma)) <= 0.0) {
						int flag_way_to_board_each = 1;
						if (border_c == 2) {
							flag_way_to_board_each = pores[j].road_to_board(pores);
							for (auto& pore : pores)
							{
								pore.set_path_to_border(0);
							}
						}
						if (flag_way_to_board_each) {
							pores_energy.push_back(energy_out);
							pores_energy_pos.push_back(j);
						}
					}
				}
			}

			for (int& j : zap_num)
			{
				for (auto poreId : pores[j].get_v_neighbor_num()) {
					if (pores[poreId].get_filled() && !pores[poreId].get_border()) {
						std::vector<int> tmp_vec = {j, poreId };
						std::ranges::sort(tmp_vec);
						if (!find_bool(pores_energy_pos_multi, tmp_vec)) {
							int flag_way_to_board_pressure = 1;
							if (border_c == 1) {
								flag_way_to_board_pressure = pores[j].get_way_to_board();
							}
							if (flag_way_to_board_pressure) {
								if ((energy_out = rabotas_out_multi(pores, tmp_vec, pres, d_sigma, sigma)) <= 0.0) {
									int flag_way_to_board_each = 1;
									if (border_c == 2) {
										flag_way_to_board_each = pores[j].road_to_board(pores);
										for (auto& pore : pores)
										{
											pore.set_path_to_border(0);
										}
									}
									if (flag_way_to_board_each) {
										pores_energy_multi.push_back(energy_out);
										pores_energy_pos_multi.push_back(tmp_vec);
									}
								}
							}
						}
					}
				}
			}

			//pores_Energy_multi = {};
			if (!pores_energy.empty() || !pores_energy_multi.empty()) {
				switch (outflow_m) {
				case 0:
				{
					for (int pores_energy_po : pores_energy_pos)
					{
						pos_e_p = pores_energy_po;
						pores[pos_e_p].set_emptied(1);
						pores[pos_e_p].set_way_to_board(0);
						pores[pos_e_p].set_filled(0);
						volume_filled -= ((4.0 * M_PI * pow(pores[pos_e_p].get_rad(), 3.0)) / 3.0);
						auto it = std::ranges::find(zap_num, pos_e_p);
						int pos_zap_num;
						pos_zap_num = std::distance(begin(zap_num), it);
						zap_num.erase(zap_num.begin() + pos_zap_num);
						out_num_new++;
						porous_neighbours.push_back(pores[pos_e_p].filled_now(pores));
					}

					for (const auto& energy_pos : pores_energy_pos_multi) {
						for (int pores_energy_po : energy_pos)
						{
							pos_e_p = pores_energy_po;
							pores[pos_e_p].set_emptied(1);
							pores[pos_e_p].set_way_to_board(0);
							pores[pos_e_p].set_filled(0);
							if (auto it = std::ranges::find(zap_num, pos_e_p); it != std::end(zap_num)) {
								volume_filled -= ((4.0 * M_PI * pow(pores[pos_e_p].get_rad(), 3.0)) / 3.0);
								int pos_zap_num;
								pos_zap_num = std::distance(begin(zap_num), it);
								zap_num.erase(zap_num.begin() + pos_zap_num);
								out_num_new++;
								porous_neighbours.push_back(pores[pos_e_p].filled_now(pores));
							}
						}

					}
					break;
				}
				case 1:
				{
					std::random_device rd;
					std::mt19937 mt(rd());
					std::uniform_int_distribution<int> dist(0, pores_energy.size() + pores_energy_multi.size() - 1);
					pos_e = dist(mt);
					if (pos_e < pores_energy.size()) {
						pos_e_p = pores_energy_pos[pos_e];
						pores[pos_e_p].set_emptied(1);
						pores[pos_e_p].set_way_to_board(0);
						pores[pos_e_p].set_filled(0);
						volume_filled -= ((4.0 * M_PI * pow(pores[pos_e_p].get_rad(), 3.0)) / 3.0);
						auto it = std::ranges::find(zap_num, pos_e_p);
						int pos_zap_num = std::distance(begin(zap_num), it);
						zap_num.erase(zap_num.begin() + pos_zap_num);
						out_num_new++;
						porous_neighbours.push_back(pores[pos_e_p].filled_now(pores));
					}
					else {
						pos_e -= pores_energy.size();
						for (auto e_p : pores_energy_pos_multi[pos_e]) {
							pores[e_p].set_emptied(1);
							pores[e_p].set_way_to_board(0);
							pores[e_p].set_filled(0);
							volume_filled -= ((4.0 * M_PI * pow(pores[e_p].get_rad(), 3.0)) / 3.0);
							auto it = std::ranges::find(zap_num, e_p);
							int pos_zap_num;
							pos_zap_num = std::distance(begin(zap_num), it);
							zap_num.erase(zap_num.begin() + pos_zap_num);
							out_num_new++;
							porous_neighbours.push_back(pores[e_p].filled_now(pores));
						}
					}
					//std::vector<double>::iterator result = select_randomly(pores_Energy.begin(), pores_Energy.end());
					//pos_e = std::distance(begin(pores_Energy), result);
					

					break;
				}
				case 2:
				{
					if (!pores_energy.empty()) {
						std::vector<double>::iterator result;
						result = std::ranges::min_element(pores_energy);
						pos_e = std::distance(begin(pores_energy), result);
						e_level = pores_energy[pos_e];
					}
					else {
						e_level = 100;
					}
					if (!pores_energy_multi.empty()) {
						std::vector<double>::iterator result_multi;
						result_multi = std::ranges::min_element(pores_energy_multi);
						pos_e_multi = std::distance(begin(pores_energy_multi), result_multi);
						e_level_multi = pores_energy_multi[pos_e_multi];
					}
					else {
						e_level_multi = 100;
					}

					if (e_level < e_level_multi) {
						pos_e_p = pores_energy_pos[pos_e];
						pores[pos_e_p].set_emptied(1);
						pores[pos_e_p].set_way_to_board(0);
						pores[pos_e_p].set_filled(0);
						volume_filled -= ((4.0 * M_PI * pow(pores[pos_e_p].get_rad(), 3.0)) / 3.0);
						auto it = std::ranges::find(zap_num, pos_e_p);
						int pos_zap_num;
						pos_zap_num = std::distance(zap_num.begin(), it);
						zap_num.erase(zap_num.begin() + pos_zap_num);
						out_num_new++;
						porous_neighbours.push_back(pores[pos_e_p].filled_now(pores));
					}
					else {
						//std::cout << "e_level_multi: "<< e_level_multi<<" e_level: "<< e_level << " pres: " << pres << " zap_num: "<<zap_num.size() << std::endl;
						for (auto e_p : pores_energy_pos_multi[pos_e_multi]) {
							//summmm++;
							pores[e_p].set_emptied(1);
							pores[e_p].set_way_to_board(0);
							pores[e_p].set_filled(0);
							volume_filled -= ((4.0 * M_PI * pow(pores[e_p].get_rad(), 3.0)) / 3.0);
							auto it = std::ranges::find(zap_num, e_p);
							int pos_zap_num;
							pos_zap_num = std::distance(zap_num.begin(), it);
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
							porous_neighbours.push_back(pores[e_p].filled_now(pores));
						}
					}
					break;
				}
				default:
					break;
				}
				if (border_c == 2) {
					for (int i : zap_num)
					{
						pores[i].set_way_to_board(0);
					}
				}
				/*
				//std::vector<double>::iterator result = select_randomly(pores_Energy.begin(), pores_Energy.end());
				std::vector<double>::iterator result = std::min_element(begin(pores_Energy), end(pores_Energy));

				pos_e = std::distance(begin(pores_Energy), result);
				pos_e_p = pores_Energy_pos[pos_e];
				pores[pos_e_p].set_emptied(1);
				pores[pos_e_p].set_filled(0);
				volume_filled -= ((4.0 * M_PI * pow(pores[pos_e_p].get_Rad(), 3.0)) / 3.0);
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
		for (int i : zap_num)
		{
			pores[i].set_way_to_board(0);
		}
		volume_filled_graph_out.push_back(volume_filled);
		pressure_graph_out.push_back(pres / 101500.0);
	}
	double max_volume_filled;
	max_volume_filled = *std::ranges::max_element(volume_filled_graph_out);
	for (double& i : volume_filled_graph_out)
	{
		i = i / max_volume_filled;
	}
	

	//std::cout << search_time << "\n";
	/*
	std::cout << porous_neighbours.size() << std::endl;
	for (auto object : porous_neighbours)
	{
		std::cout << object << " ";
	}
	std::cout << std::endl;
	*/
}

void generate_body(const int n, const std::vector<int> core_v, const std::vector<int> d_sigma_v, const int i_min, const int i_max, const std::string connection_string, const int body_id, std::vector<int> percent_for_del_vec) {
	mtx.lock();
	std::cout << "Generate id: " << std::this_thread::get_id() << " started" << std::endl;
	mtx.unlock();
	for (int i = i_min; i < i_max; i++) {
		for (const int j : d_sigma_v)
		{
			for (const auto percent_for_del : percent_for_del_vec) {
				std::string table_name = "body_core_" + std::to_string(core_v[i]) + "_DS_" + std::to_string(j) + "_ID_" + std::to_string(body_id)+ "_PFD_"+ std::to_string(percent_for_del);
				std::string graph_id = std::to_string(core_v[i]) + std::to_string(j) + std::to_string(body_id) + std::to_string(percent_for_del);
				std::vector<pore> pores;
				std::vector<double> volume_filled_graph;
				std::vector<double> pressure_graph;
				double volume_filled = 0;
				std::vector<int> porous_neighbours;
				filling(pores, volume_filled_graph, pressure_graph, volume_filled, n, core_v[i], j, percent_for_del, porous_neighbours);
				std::string connection_string_graphs = "host=localhost port=1768 dbname=graphs user=postgres password =adu202121";
				insert_into_table_graph(graph_id, volume_filled_graph, pressure_graph, "graph_in", connection_string_graphs);
				create_table(table_name, connection_string);
				insert_into_table(pores, table_name, connection_string);
				connection_string_graphs = "host=localhost port=1768 dbname=graphs user=postgres password =adu202121";
				insert_into_table_graph_li(graph_id, porous_neighbours, connection_string_graphs);
			}
		}
	}
	mtx.lock();
	std::cout << "Generate id: " << std::this_thread::get_id() << " finished" << std::endl;
	mtx.unlock();
}

void empty_body(std::vector<int> core_v, std::vector<int> d_sigma_v, int i_min, int i_max, std::string connectionString, std::vector<int> body_id, int border_C, int outflow_M, std::vector<int> percent_for_del_vec) {
	mtx.lock();
	std::cout << "Empty id: " << std::this_thread::get_id() << " " << border_C << " " << outflow_M << " started" << std::endl;
	mtx.unlock();
	for (const int k : body_id)
	{
		for (int i = i_min; i < i_max; i++) {
			for (const int j : d_sigma_v)
			{
				for (const auto percent_for_del : percent_for_del_vec) {
					std::string table_name = "body_core_" + std::to_string(core_v[i]) + "_DS_" + std::to_string(j) + "_ID_" + std::to_string(
						k) + "_PFD_" + std::to_string(percent_for_del);
					std::string graph_id = std::to_string(core_v[i]) + std::to_string(j) + std::to_string(k) + std::to_string(border_C) + std::to_string(outflow_M) + std::to_string(percent_for_del);
					std::vector<pore> pores;
					select_from_table(pores, table_name, connectionString);
					double volume_filled = 0;

					for (auto& pore : pores)
					{
						if (pore.get_filled() && !(pore.get_border())) {
							volume_filled += ((4.0 * M_PI * pow(pore.get_rad(), 3.0)) / 3.0);
						}
					}
					std::vector<double> volume_filled_graph_out;
					std::vector<double> pressure_graph_out;
					std::vector<int> porous_neighbours;
					leakage(pores, volume_filled_graph_out, pressure_graph_out, volume_filled, core_v[i], j, border_C, outflow_M, porous_neighbours);

					std::string connection_string_graphs = "host=localhost port=1768 dbname=graphs user=postgres password =adu202121";
					insert_into_table_graph(graph_id, volume_filled_graph_out, pressure_graph_out, "graph_out", connection_string_graphs);
					connection_string_graphs = "host=localhost port=1768 dbname=graphs user=postgres password =adu202121";
					insert_into_table_graph_li(graph_id, porous_neighbours, connection_string_graphs);
				}
			}
		}
	}
	mtx.lock();
	std::cout << "Empty id: " << std::this_thread::get_id() << " " << border_C << " " << outflow_M << " finished" << std::endl;
	mtx.unlock();
}


int main() {
	///
	///
	/// ѕовысить core что позволит уменьшить количство сеседей, но сильно уменьшить d_sigma
	/// 
	/// ѕосле генерации тела по гаусовому закону, пройти по порам и у пор с большим количеством соседей заменить радиус случайной величининой из  эксперимента
	/// 
	///
	///
	/// 
	/// 
	///
	simple_timer st;
	SetConsoleOutputCP(1251);
	std::string connectionString = "host=localhost port=1768 dbname=porous_body user=postgres password =adu202121";

	//pqxx::result response;
	/*
	pqxx::result response = worker.exec("SELECT * FROM body_1 WHERE id=3");

	for (int i = 0; i < response.size(); i++)
	{
		for (int j = 0; j < response[i].size(); j++)
		{
			str += to_string(response[i][j])+" ";
		}

		std::cout <<str<< std::endl;

	}
*/

// border_C: 0 - не учитываем, 1 - расчЄт после измениени€ давлени€, 2 - расчЄт дл€ каждой
// outflow_M: 0 - все, 1 - случайна€, 2 - с минимальной энерги€ей
//std::vector<int> body_id = { 1,2,3,4,5,6,7,8,9,10 };
	std::vector<int> body_id = { 80 };
	std::vector<int> border_c = { 0,1,2 };
	std::vector<int> outflow_m = { 0,1,2 };
	std::vector<int> core_v;
	std::vector<int> d_sigma_v;
	for (int i = 0; i <= 10; i++) {
		core_v.push_back(93 + 2 * i);
	}
	for (int i = 0; i <= 0; i++) {
		d_sigma_v.push_back(9 + 2 * i);
	}

	std::vector<int> percent_for_del_vec = { 10 };

	std::cout << "ѕор на границе: ";
	int n;
	std::cin >> n;

	/////////////////////////////////////////////////
	std::thread th0_in(generate_body, n, core_v, d_sigma_v, 0 * core_v.size() / 11, 1 * core_v.size() / 11, connectionString, body_id[0], percent_for_del_vec);
	
	//std::thread th1_in(generate_body, N, core_v, d_sigma_v, 1 * core_v.size() / 11, 2 * core_v.size() / 11, connectionString, body_id[0], percent_for_del_vec);
	//std::thread th2_in(generate_body, N, core_v, d_sigma_v, 2 * core_v.size() / 11, 3 * core_v.size() / 11, connectionString, body_id[0], percent_for_del_vec);
	/*
	std::thread th3_in(generate_body, N, core_v, d_sigma_v, 3 * core_v.size() / 11, 4 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th4_in(generate_body, N, core_v, d_sigma_v, 4 * core_v.size() / 11, 5 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th5_in(generate_body, N, core_v, d_sigma_v, 5 * core_v.size() / 11, 6 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th6_in(generate_body, N, core_v, d_sigma_v, 6 * core_v.size() / 11, 7 * core_v.size() / 11, connectionString, body_id[0]);
	*/
	/*
	std::thread th7_in(generate_body, N, core_v, d_sigma_v, 7 * core_v.size() / 11, 8 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th8_in(generate_body, N, core_v, d_sigma_v, 8 * core_v.size() / 11, 9 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th9_in(generate_body, N, core_v, d_sigma_v, 9 * core_v.size() / 11, 10 * core_v.size() / 11, connectionString, body_id[0]);
	std::thread th10_in(generate_body, N, core_v, d_sigma_v, 10 * core_v.size() / 11, 11 * core_v.size() / 11, connectionString, body_id[0]);
	*/

	/////////////////////////////////////////////////
	th0_in.join();
	
	//th1_in.join();
	//th2_in.join();
	/*
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
	
	std::thread th00_out(empty_body, core_v, d_sigma_v, 0 * core_v.size() / 11, 1 * core_v.size() / 11, connectionString, body_id, border_c[2], outflow_m[0], percent_for_del_vec);
	
	//std::thread th01_out(empty_body, core_v, d_sigma_v, 1 * core_v.size() / 11, 2 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0], percent_for_del_vec);
	//std::thread th02_out(empty_body, core_v, d_sigma_v, 2 * core_v.size() / 11, 3 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0], percent_for_del_vec);
	/*
	std::thread th03_out(empty_body, core_v, d_sigma_v, 3 * core_v.size() / 11, 4 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th04_out(empty_body, core_v, d_sigma_v, 4 * core_v.size() / 11, 5 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th05_out(empty_body, core_v, d_sigma_v, 5 * core_v.size() / 11, 6 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th06_out(empty_body, core_v, d_sigma_v, 6 * core_v.size() / 11, 7 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	*/
	/*
	std::thread th07_out(empty_body, core_v, d_sigma_v, 7 * core_v.size() / 11, 8 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th08_out(empty_body, core_v, d_sigma_v, 8 * core_v.size() / 11, 9 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th09_out(empty_body, core_v, d_sigma_v, 9 * core_v.size() / 11, 10 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	std::thread th010_out(empty_body, core_v, d_sigma_v, 10 * core_v.size() / 11, 11 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[0]);
	*/
	
	//////////////////////////////////////////
	
	th00_out.join();
	
	//th01_out.join();
	//th02_out.join();
	/*
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
	
	//std::thread th10_out(empty_body, core_v, d_sigma_v, 0 * core_v.size() / 11, 1 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1], percent_for_del_vec);
	/*
	std::thread th11_out(empty_body, core_v, d_sigma_v, 1 * core_v.size() / 11, 2 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th12_out(empty_body, core_v, d_sigma_v, 2 * core_v.size() / 11, 3 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th13_out(empty_body, core_v, d_sigma_v, 3 * core_v.size() / 11, 4 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th14_out(empty_body, core_v, d_sigma_v, 4 * core_v.size() / 11, 5 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th15_out(empty_body, core_v, d_sigma_v, 5 * core_v.size() / 11, 6 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th16_out(empty_body, core_v, d_sigma_v, 6 * core_v.size() / 11, 7 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	*/
	/*
	std::thread th17_out(empty_body, core_v, d_sigma_v, 7 * core_v.size() / 11, 8 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th18_out(empty_body, core_v, d_sigma_v, 8 * core_v.size() / 11, 9 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th19_out(empty_body, core_v, d_sigma_v, 9 * core_v.size() / 11, 10 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
	std::thread th110_out(empty_body, core_v, d_sigma_v, 10 * core_v.size() / 11, 11 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[1]);
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
	
	//std::thread th20_out(empty_body, core_v, d_sigma_v, 0 * core_v.size() / 11, 1 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2], percent_for_del_vec);
	/*
	std::thread th21_out(empty_body, core_v, d_sigma_v, 1 * core_v.size() / 11, 2 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th22_out(empty_body, core_v, d_sigma_v, 2 * core_v.size() / 11, 3 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th23_out(empty_body, core_v, d_sigma_v, 3 * core_v.size() / 11, 4 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th24_out(empty_body, core_v, d_sigma_v, 4 * core_v.size() / 11, 5 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th25_out(empty_body, core_v, d_sigma_v, 5 * core_v.size() / 11, 6 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th26_out(empty_body, core_v, d_sigma_v, 6 * core_v.size() / 11, 7 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	*/
	/*
	std::thread th27_out(empty_body, core_v, d_sigma_v, 7 * core_v.size() / 11, 8 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th28_out(empty_body, core_v, d_sigma_v, 8 * core_v.size() / 11, 9 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th29_out(empty_body, core_v, d_sigma_v, 9 * core_v.size() / 11, 10 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	std::thread th210_out(empty_body, core_v, d_sigma_v, 10 * core_v.size() / 11, 11 * core_v.size() / 11, connectionString, body_id, border_C[2], outflow_M[2]);
	*/

	//////////////////////////////////////////
	//th20_out.join();
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
