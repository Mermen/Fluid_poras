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
#include  <iterator>
#include <pqxx/pqxx>
#include <cstdio>
#include <windows.h>
//using namespace std;

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

double linear(std::vector<double>& x, std::vector<double>& y, std::vector<double>& k, std::vector<double>& b, double X, int n) {
    // за пределами точек
    if (X <= x[0])
        return y[0];
    else if (X >= x[n - 1])
        return y[n - 1];

    // между точками
    /*
    for (int i = 0; i < n - 1; ++i)
        if (X > x[i] && X < x[i + 1])
            return k[i] * X + b[i];
    */
    int middle = 0;
    int left = 0;
    int right = n;

    while (1) {
        middle = (left + right) / 2;
        if (X < x[middle]) {
            right = middle;
        }
        else if (X > x[middle]) {
            left = middle;
        }
        if ((left + 1) == right) {
            return k[left] * X + b[left];
        }

    }


    // в точках
    for (int i = 0; i < n - 1; ++i)
        if (X == x[i])
            return y[i];


    return -1;	// ошибка
}

class Pora {
    double X;
    double Y;
    double Z;
    double Rad;
    int border;
    int filled;
    int emptied;
    int path_to_border;
    int way_to_board;
    std::vector<int> neighbors_num;
    std::vector<double> neighbors_area;
public: Pora(double x, double y, double z, double Rad, int border, int filled, int emptied, int path_to_border, int way_to_board) {
    this->X = x;
    this->Y = y;
    this->Z = z;
    this->Rad = Rad;
    this->border = border;
    this->filled = filled;
    this->emptied = emptied;
    this->path_to_border = path_to_border;
    this->way_to_board = way_to_board;
}
      
      double get_X() {
          return this->X;
      }
      double get_Y() {
          return this->Y;
      }
      double get_Z() {
          return this->Z;
      }
      double get_Rad() {
          return this->Rad;
      }
      int get_border() {
          return this->border;
      }
      int get_filled() {
          return this->filled;
      }
      int get_emptied() {
          return this->emptied;
      }
      int get_path_to_border() {
          return this->path_to_border;
      }
      int get_way_to_board() {
          return this->way_to_board;
      }
      int get_neighbor_Num(int i) {
          return this->neighbors_num[i];
      }
      std::vector<int> get_V_neighbor_Num() {
          return this->neighbors_num;
      }
      double get_neighbor_Area(int i) {
          return this->neighbors_area[i];
      }
      std::vector<double> get_V_neighbor_Area() {
          return this->neighbors_area;
      }

      void set_border(int border) {
          this->border = border;
      }
      void set_filled(int filled) {
          this->filled = filled;
      }
      void set_emptied(int emptied) {
          this->emptied = emptied;
      }
      void set_neighbor_Num(int k) {
          this->neighbors_num.push_back(k);
      }
      void set_neighbor_Area(double area) {
          this->neighbors_area.push_back(area);
      }
      void set_path_to_border(int path_to_border) {
          this->path_to_border = path_to_border;
      }
      void set_way_to_board(int way_to_board) {
          this->way_to_board = way_to_board;
      }
      void set_neighbors_num(std::vector<int> neighbors_num) {
          this->neighbors_num = neighbors_num;
      }
      void set_neighbors_area(std::vector<double> neighbors_area) {
          this->neighbors_area = neighbors_area;
      }

      int find_k(int k) {
          for (int i = 0; i < this->neighbors_num.size(); i++)
          {
              if (this->neighbors_num[i] == k)
              {
                  return 1;

              }
          }
          return 0;

      }

      int road_to_board(std::vector<Pora>& poras) {

          if (this->border)
          {
              return 1;
          }
          else if (this->emptied) {
              return 0;
          }
          else if (this->path_to_border)
          {
              return 0;
          }
          else
          {
              this->path_to_border = 1;
              int road = 0;
              for (int i = 0; i < this->neighbors_num.size(); i++)
              {
                  road += poras[neighbors_num[i]].road_to_board(poras);
              }
              return (road > 0);
          }
      }
};

std::string to_string_pora(Pora pora) {
    std::string str;
    std::string str_num="";
    std::string str_area="";
    std::ostringstream ostr_num;
    std::ostringstream ostr_area;
    std::vector<int> vec_num;
    std::vector<double> vec_area;
    vec_num = pora.get_V_neighbor_Num();
    vec_area = pora.get_V_neighbor_Area();

    str_num += '\'';
    str_num += '{';
    if (!vec_num.empty())
    {
        // Convert all but the last element to avoid a trailing ","
        std::copy(vec_num.begin(), vec_num.end() - 1,
            std::ostream_iterator<int>(ostr_num, ","));

        // Now add the last element with no delimiter
        ostr_num << vec_num.back();
    }
    str_num += ostr_num.str();
    str_num += '}';
    str_num += '\'';

    str_area += '\'';
    str_area += '{';
    if (!vec_area.empty())
    {
        // Convert all but the last element to avoid a trailing ","
        std::copy(vec_area.begin(), vec_area.end() - 1,
            std::ostream_iterator<int>(ostr_area, ","));

        // Now add the last element with no delimiter
        ostr_area << vec_area.back();
    }
    str_area += ostr_area.str();
    str_area += '}';
    str_area += '\'';

    str = std::to_string(pora.get_X()) + "," + std::to_string(pora.get_Y()) + "," + std::to_string(pora.get_Z()) + "," + std::to_string(pora.get_Rad()) + ","
        + std::to_string(pora.get_border()) + "," + std::to_string(pora.get_filled()) + "," + std::to_string(pora.get_emptied()) + ","
        + std::to_string(pora.get_path_to_border()) + "," + std::to_string(pora.get_way_to_board()) + ","
        + str_num + "," + str_area;
    return str;
}

void insert_into_table(std::vector<Pora>& poras, std::string table_name, std::string connectionString) {
    pqxx::connection connectionObject(connectionString.c_str());
    pqxx::work worker(connectionObject);
    worker.exec("set client_encoding='win1251'");
    for (size_t i = 0; i < poras.size(); i++)
    {
        worker.exec("INSERT INTO " + table_name + " values("+std::to_string(i)+","
        + to_string_pora(poras[i])+");");
    }
    worker.commit();
    connectionObject.close();
}

void select_from_table(std::vector<Pora>& poras, std::string table_name, std::string connectionString) {
    std::string str = "";
    std::string tmp;
    pqxx::connection connectionObject(connectionString.c_str());
    pqxx::work worker(connectionObject);
    worker.exec("set client_encoding='win1251'");
    pqxx::result response = worker.exec("SELECT * FROM "+table_name);

    for (size_t i = 0; i < response.size(); i++)
    {
        poras.push_back(Pora(std::stod(to_string(response[i][1])), std::stod(to_string(response[i][2])), std::stod(to_string(response[i][3])), std::stod(to_string(response[i][4])),
            std::stoi(to_string(response[i][5])), std::stoi(to_string(response[i][6])), std::stoi(to_string(response[i][7])), std::stoi(to_string(response[i][8])), std::stoi(to_string(response[i][9]))));
        std::stringstream ss0(to_string(response[i][10]).substr(1, to_string(response[i][10]).size() - 2));
        std::vector<int> neighbors_num;
        while (std::getline(ss0, tmp, ',')) {
            neighbors_num.push_back(std::stoi(tmp));
        }
        poras[i].set_neighbors_num(neighbors_num);

        std::stringstream ss1(to_string(response[i][11]).substr(1, to_string(response[i][11]).size() - 2));
        std::vector<double> neighbors_area;
        while (std::getline(ss1, tmp, ',')) {
            neighbors_area.push_back(std::stod(tmp));
        }
        poras[i].set_neighbors_area(neighbors_area);

        /*
        str += to_string(response[i][0]) + " ";
        for (size_t j = 1; j < 5; j++)
        {
            str += to_string(response[i][j]) + " ";
        }
        for (size_t j = 5; j < 10; j++)
        {
            str += to_string(response[i][j]) + " ";
        }

        std::stringstream ss0(to_string(response[i][10]).substr(1, to_string(response[i][10]).size() - 2));
        std::vector<int> neighbors_num;
        while (std::getline(ss0, tmp, ',')) {
            neighbors_num.push_back(std::stoi(tmp));
        }
        for (size_t j = 0; j < neighbors_num.size(); j++)
        {
            str += std::to_string(neighbors_num[j])+" ";
        }

        std::stringstream ss1(to_string(response[i][11]).substr(1, to_string(response[i][11]).size() - 2));
        std::vector<double> neighbors_area;
        while (std::getline(ss1, tmp, ',')) {
            neighbors_area.push_back(std::stod(tmp));
        }
        for (size_t j = 0; j < neighbors_area.size(); j++)
        {
            str += std::to_string(neighbors_area[j]) + " ";
        }
        std::cout << str << std::endl;
        str = "";
        */
    }
    worker.commit();
    connectionObject.close();
}

void create_table(std::string table_name, std::string connectionString) {
    std::string str = "";
    std::string tmp;
    pqxx::connection connectionObject(connectionString.c_str());
    pqxx::work worker(connectionObject);
    worker.exec("set client_encoding='win1251'");
    pqxx::result response = worker.exec("CREATE TABLE " + table_name
    +" (id integer, x double precision, y double precision, z double precision, rad double precision, border integer, filled integer, emptied integer, path_to_border integer, way_to_border integer, neighbors_num integer[], neighbors_area double precision[]);");
    worker.commit();
    connectionObject.close();
}

void Generation_Poras_gauss(double Rad, double core, std::vector<Pora>& poras, double N, double sii) {
    std::random_device rd_R;
    std::normal_distribution<> uid_R(Rad, sii);
    bool kj = 0;
    double Rad_tmp = 0;
    double pot_por_rad = 0;
    int start_time = clock();
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                poras.push_back(Pora(Rad * core / 2.0 + 1.0 * i * (Rad + Rad) * core, Rad * core / 2.0 + 1.0 * j * (Rad + Rad) * core, Rad * core / 2.0 + 1.0 * k * (Rad + Rad) * core, uid_R(rd_R), 0, 0, 0, 0, 0));
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
                poras.push_back(Pora(Rad * core / 2.0 + 1.0 * i * (Rad + Rad) * core, Rad * core / 2.0 + 1.0 * j * (Rad + Rad) * core, Rad * core / 2.0 + 1.0 * k * (Rad + Rad) * core, Rad_tmp, 0, 0, 0, 0, 0));
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
double rabotas_out(std::vector<Pora>& poras, int i, double pres, double DSIGMA, double SIGMA) {
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
    work = pres * 4.0 * pow(10.0, -9.0) * M_PI / 3.0 * pow(poras[i].get_Rad(), 3.0)
        - (4.0 * M_PI * poras[i].get_Rad() * poras[i].get_Rad() - SM) * DSIGMA
        - SIGMA * (surf_empty - SMZ);
    return work;
}

int get_neibor_num(std::vector<Pora>& poras, int i) {
    int num = 0;
    for (int j = 0; j < poras[i].get_V_neighbor_Num().size(); j++)
    {
        if (poras[poras[i].get_neighbor_Num(j)].get_filled())
        {
            num++;
        }
    }
    return num;
}


void filling(std::vector<Pora>& poras, std::vector<double>& volume_filled_graph, std::vector<double>& pressure_graph, double& volume_filled_sum, int N, double core, double DSIGMA, std::string table_name) {

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

    Generation_Poras_experiment(x, y, interp_k, interp_b, Rad, core, poras, N);
    //Generation_Poras_gauss(Rad, core, poras, N, sii);


    //int aa;
    //cin >> aa;
    std::vector<int> bound_1;
    std::vector<int> bound_2;
    std::vector<int> bound_3;
    std::vector<int> bound_4;
    std::vector<int> bound_5;
    std::vector<int> bound_6;
    std::vector<int> bound;
    for (int i = 0; i < poras.size(); i++) {
        if (poras[i].get_X() == Rad * core / 2.0) {
            poras[i].set_border(1);
            bound_1.push_back(i);
        }
        else if (poras[i].get_X() == Rad * core / 2.0 + (N - 1) * (Rad + Rad) * core) {
            poras[i].set_border(1);
            bound_2.push_back(i);
        }
        if (poras[i].get_Y() == Rad * core / 2.0) {
            poras[i].set_border(1);
            bound_3.push_back(i);
        }
        else if (poras[i].get_Y() == Rad * core / 2.0 + (N - 1) * (Rad + Rad) * core) {
            poras[i].set_border(1);
            bound_4.push_back(i);
        }
        if (poras[i].get_Z() == Rad * core / 2.0) {
            poras[i].set_border(1);
            bound_5.push_back(i);
        }
        else if (poras[i].get_Z() == Rad * core / 2.0 + (N - 1) * (Rad + Rad) * core) {
            poras[i].set_border(1);
            bound_6.push_back(i);
        }
        if (poras[i].get_border()) {
            bound.push_back(i);
        }

    }
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
            else if (k <= (poras.size() - 1) && k >= 0)
            {
                kappa = pow((poras[i].get_X() - poras[k].get_X()) * (poras[i].get_X() - poras[k].get_X()) +
                    (poras[i].get_Y() - poras[k].get_Y()) * (poras[i].get_Y() - poras[k].get_Y()) +
                    (poras[i].get_Z() - poras[k].get_Z()) * (poras[i].get_Z() - poras[k].get_Z()), 0.5);
                if (kappa <= (poras[i].get_Rad() + poras[k].get_Rad()))
                {
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

    for (int i = 0; i < poras.size(); i++)
    {
        if (poras[i].get_border())
        {
            poras[i].set_filled(1);
            zap_num.push_back(i);
        }
    }
    int zap_num_new = 0;
    double volume_filled = 0;
    int time_zap_start = clock();
    for (double pres = pressure_start; pres < pressure_end; pres += pressure_step) {
        zap_num_new = 1;
        std::cout << pres << "\n";
        while (zap_num_new) {
            zap_num_new = 0;

            for (int j = 0; j < zap_num.size(); j++) {
                for (int k = 0; k < poras[zap_num[j]].get_V_neighbor_Num().size(); k++) {
                    if (!(poras[poras[zap_num[j]].get_neighbor_Num(k)].get_filled())) {
                        if (rabotas_in(poras, poras[zap_num[j]].get_neighbor_Num(k), pres, DSIGMA, SIGMA))
                        {
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
    std::cout << N << "\n";
    std::cout << time_zap_serch << "\n";
    double max_volume_filled;
    max_volume_filled = *max_element(volume_filled_graph.begin(), volume_filled_graph.end());
    volume_filled_sum = max_volume_filled;
    for (int i = 0; i < volume_filled_graph.size(); i++)
    {
        volume_filled_graph[i] = volume_filled_graph[i] / max_volume_filled;
    }


    int end_time = clock();
    int search_time = end_time - start_time;
    std::cout << search_time << "\n";
    /*
    ofstream fout("log.txt");
    fout << "time_zap" << "		" << time_zap_serch << "\n";
    fout << "time_zap_work" << "		" << search_time << "\n";
    for (int i = 0; i < poras.size(); i++) {
        fout << poras[i].get_Rad();
        if (!(i == (poras.size() - 1)))
            fout << "\n";
    }
    ofstream bout("bound.txt");
    for (int i = 0; i < bound.size(); i++) {
        bout << bound[i];
        if (!(i == (bound.size() - 1)))
            bout << "\n";
    }

    ofstream bout_area("bound_area.txt");
    for (int i = 0; i < poras.size(); i++)
    {
        bout_area << i << "	-	" << poras[i].get_Rad() << "\n";

        for (int j = 0; j < poras[i].get_V_neighbor_Num().size(); j++)
        {
            bout_area << "	" << poras[i].get_neighbor_Num(j) << " - " << poras[i].get_neighbor_Area(j);
            if (!((i == (poras.size() - 1 ))&&(j == (poras[i].get_V_neighbor_Num().size() - 1))))
                bout_area << "\n";
        }
        if (!(i == (poras.size() - 1)))
            bout_area << "\n";
    }

    bout_area.close();
    fin.close();
    fout.close();
    bout.close();
    */

    std::ofstream gout("graph_in_"+ table_name +".csv");
    gout << std::setprecision(16) << std::fixed;
    for (int i = 0; i < volume_filled_graph.size(); i++) {
        gout << std::setprecision(16) << volume_filled_graph[i] << "," << std::setprecision(16) << pressure_graph[i];
        if (!(i == (volume_filled_graph.size() - 1)))
            gout << "\n";
    }
    gout.close();




}

void leakage(std::vector<Pora>& poras, std::vector<double>& volume_filled_graph_out, std::vector<double>& pressure_graph_out, double volume_filled, double core, double DSIGMA, std::string table_name) {
    int start_time = clock();
    std::vector<int> zap_num;
    for (int i = 0; i < poras.size(); i++)
    {
        if (poras[i].get_filled())
        {
            zap_num.push_back(i);
        }
    }
    double pressure_start = 0;
    double pressure_end = 34000001;
    double pressure_step = 100000;
    //double DSIGMA = 40.0 / 1000.0; //DSIGMA поверхностна€ энерги€
    double SIGMA = 75.3 / 1000.0;//SIGMA поверхностна€ энерги€ жидкость-газ
    double out_num_new = 0;
    std::vector<double> poras_Energy(0);
    std::vector<int> poras_Energy_pos(0);
    double energy_out = 0.0;
    int pos_e = 0;
    int pos_e_p = 0;
    for (double pres = pressure_end - 1; pres > pressure_start - 1; pres -= pressure_step) {
        out_num_new = 1;
        std::cout << pres << "\n";
        /*
        for (int i = 0; i < zap_num.size(); i++)
        {
            if (poras[zap_num[i]].road_to_board(poras)) {
                for (int i = 0; i < poras.size(); i++)
                {
                    poras[i].set_path_to_border(0);
                }
                poras[zap_num[i]].set_way_to_board(1);
            }
            else
            {
                for (int i = 0; i < poras.size(); i++)
                {
                    poras[i].set_path_to_border(0);
                }
            }
        }
        */
        while (out_num_new) {
            out_num_new = 0;
            //for (int i = 1; i <= 26; i++) {
            poras_Energy.clear();
            poras_Energy.shrink_to_fit();
            poras_Energy_pos.clear();
            poras_Energy_pos.shrink_to_fit();
            for (int j = 0; j < zap_num.size(); j++) {
                if (!(poras[zap_num[j]].get_border())) {
                    if (!(poras[zap_num[j]].get_emptied())) {
                        //if (poras[zap_num[j]].get_way_to_board()) {
                        if ((energy_out = rabotas_out(poras, zap_num[j], pres, DSIGMA, SIGMA)) <= 0.0) {
                            
                            if (poras[zap_num[j]].road_to_board(poras))
                            {
                                for (int i = 0; i < poras.size(); i++)
                                {
                                    poras[i].set_path_to_border(0);
                                }
                                
                                poras_Energy.push_back(energy_out);
                                poras_Energy_pos.push_back(zap_num[j]);
                                
                            }
                            else
                            {
                                for (int i = 0; i < poras.size(); i++)
                                {
                                    poras[i].set_path_to_border(0);
                                }
                            }
                            
                            /*
                            poras[zap_num[j]].set_emptied(1);
                            poras[zap_num[j]].set_filled(0);
                            volume_filled -= ((4.0 * M_PI * pow(poras[zap_num[j]].get_Rad(), 3.0)) / 3.0);
                            zap_num.erase(zap_num.begin() + j);
                            j = -1;
                            i = 1;
                            out_num_new++;
                            */
                        }
                       // }
                    }
                }
            }
            if (poras_Energy.size())
            {
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
                
                /*
                                    for (int k = 0; k < poras_Energy_pos.size(); k++)
                                    {
                                        pos_e_p = poras_Energy_pos[k];
                                        poras[pos_e_p].set_emptied(1);
                                        poras[pos_e_p].set_filled(0);
                                        volume_filled -= ((4.0 * M_PI * pow(poras[pos_e_p].get_Rad(), 3.0)) / 3.0);
                                        auto it = std::find(zap_num.begin(), zap_num.end(), pos_e_p);
                                        int pos_zap_num = std::distance(begin(zap_num), it);
                                        zap_num.erase(zap_num.begin() + pos_zap_num);
                                        //j = -1;
                                        //i = 0;
                                        out_num_new++;
                                    }
                */

            }
            //}

        }
        for (int i = 0; i < zap_num.size(); i++)
        {
            poras[zap_num[i]].set_way_to_board(0);
        }
        volume_filled_graph_out.push_back(volume_filled);
        pressure_graph_out.push_back(pres / 101500.0);
    }
    double max_volume_filled;
    max_volume_filled = *max_element(volume_filled_graph_out.begin(), volume_filled_graph_out.end());
    for (int i = 0; i < volume_filled_graph_out.size(); i++)
    {
        volume_filled_graph_out[i] = volume_filled_graph_out[i] / max_volume_filled;
    }
    int end_time = clock();
    int search_time = end_time - start_time;
    std::cout << search_time << "\n";
    std::ofstream gout("graph_out.csv");
    gout << std::setprecision(16) << std::fixed;
    for (int i = 0; i < volume_filled_graph_out.size(); i++) {
        gout << std::setprecision(16) << volume_filled_graph_out[i] << "," << std::setprecision(16) << pressure_graph_out[i];
        if (!(i == (volume_filled_graph_out.size() - 1)))
            gout << "\n";
    }
    gout.close();
}



int main() {
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

    std::vector<int> core_v;
    std::vector<int> DSIGMA_v;
    for (size_t i = 1; i <= 4; i++)
    {
        core_v.push_back(50 + i);
    }
    for (size_t i = 1; i <= 4; i++)
    {
        DSIGMA_v.push_back(i);
    }
    
    std::cout << "ѕор на границе: ";
    int N;
    std::cin >> N;
    for (size_t i = 0; i < core_v.size(); i++)
    {
        for (size_t j = 0; j < DSIGMA_v.size(); j++)
        {
            double core = (1.0 * core_v[i]) / 100.0;
            double DSIGMA = (1.0 * DSIGMA_v[j]) / 1000.0;
            std::string table_name = "body_core_" + std::to_string(core_v[i]) + "_DS_" + std::to_string(DSIGMA_v[j]);
            std::vector<Pora> poras;
            std::vector<double> volume_filled_graph;
            std::vector<double> pressure_graph;
            double volume_filled = 0;
            filling(poras, volume_filled_graph, pressure_graph, volume_filled, N, core, DSIGMA, table_name);
            create_table(table_name, connectionString);
            insert_into_table(poras, table_name, connectionString);
        }
    }
    std::cout << "Stop test" << std::endl;
    std::cin >> N;
    /*
    select_from_table(poras, "body_1", connectionString);
    std::cout << "select done" << std::endl;
    std::cin>>N;
    double volume_filled = 0;
    if (volume_filled == 0)
    {
        for (size_t i = 0; i < poras.size(); i++)
        {
            if (poras[i].get_filled())
            {

                volume_filled += ((4.0 * M_PI * pow(poras[i].get_Rad(), 3.0)) / 3.0);
            }
        }
    }
    std::cout << "Volume filled: " << volume_filled;
    std::cin >> N;
    std::vector<double> volume_filled_graph_out;
    std::vector<double> pressure_graph_out;
    leakage(poras, volume_filled_graph_out, pressure_graph_out, volume_filled, core, DSIGMA);
    */



}
