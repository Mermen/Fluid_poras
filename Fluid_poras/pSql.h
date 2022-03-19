#pragma once
void insert_into_table(std::vector<Pora>& poras, std::string table_name, std::string connectionString) {
    pqxx::connection connectionObject(connectionString.c_str());
    pqxx::work worker(connectionObject);
    worker.exec("set client_encoding='win1251'");
    for (size_t i = 0; i < poras.size(); i++)
    {
        worker.exec("INSERT INTO " + table_name + " values(" + std::to_string(i) + ","
            + to_string_pora(poras[i]) + ");");
    }
    worker.commit();
    connectionObject.close();
}

void insert_into_table_graph(std::string graph_id, std::vector<double> volume_filled_graph, std::vector<double> pressure_graph, std::string table_name, std::string connectionString) {
    std::string str_volume = "";
    std::string str_preassure = "";
    std::ostringstream ostr_volume;
    std::ostringstream ostr_preassure;
    std::vector<double> vec_volume;
    std::vector<double> vec_preassure;
    vec_volume = volume_filled_graph;
    vec_preassure = pressure_graph;

    str_volume += '\'';
    str_volume += '{';
    if (!vec_volume.empty())
    {
        // Convert all but the last element to avoid a trailing ","
        std::copy(vec_volume.begin(), vec_volume.end() - 1,
            std::ostream_iterator<double>(ostr_volume, ","));

        // Now add the last element with no delimiter
        ostr_volume << vec_volume.back();
    }
    str_volume += ostr_volume.str();
    str_volume += '}';
    str_volume += '\'';

    str_preassure += '\'';
    str_preassure += '{';
    if (!vec_preassure.empty())
    {
        // Convert all but the last element to avoid a trailing ","
        std::copy(vec_preassure.begin(), vec_preassure.end() - 1,
            std::ostream_iterator<double>(ostr_preassure, ","));

        // Now add the last element with no delimiter
        ostr_preassure << vec_preassure.back();
    }
    str_preassure += ostr_preassure.str();
    str_preassure += '}';
    str_preassure += '\'';


    pqxx::connection connectionObject(connectionString.c_str());
    pqxx::work worker(connectionObject);
    worker.exec("set client_encoding='win1251'");
    worker.exec("DELETE FROM " + table_name + " WHERE id=" + graph_id + ";");
    worker.exec("INSERT INTO " + table_name + " values(" + graph_id + ","
        + str_volume + "," + str_preassure + ",'url');");
    worker.commit();
    connectionObject.close();
}

void insert_into_table_graph_li(std::string graph_id, std::vector<int> porous_neibors, std::string connectionString) {
    std::string str_neibors = "";
    std::ostringstream ostr_neibors;
    std::vector<int> vec_neibors;
    vec_neibors = porous_neibors;

    str_neibors += '\'';
    str_neibors += '{';
    if (!vec_neibors.empty())
    {
        // Convert all but the last element to avoid a trailing ","
        std::copy(vec_neibors.begin(), vec_neibors.end() - 1,
            std::ostream_iterator<double>(ostr_neibors, ","));

        // Now add the last element with no delimiter
        ostr_neibors << vec_neibors.back();
    }
    str_neibors += ostr_neibors.str();
    str_neibors += '}';
    str_neibors += '\'';


    pqxx::connection connectionObject(connectionString.c_str());
    pqxx::work worker(connectionObject);
    worker.exec("set client_encoding='win1251'");
    worker.exec("DELETE FROM graph_li WHERE id=" + graph_id + ";");
    worker.exec("INSERT INTO graph_li values(" + graph_id + ","
        + str_neibors + ");");
    worker.commit();
    connectionObject.close();
}

void select_from_table(std::vector<Pora>& poras, std::string table_name, std::string connectionString) {
    std::string str = "";
    std::string tmp;
    pqxx::connection connectionObject(connectionString.c_str());
    pqxx::work worker(connectionObject);
    worker.exec("set client_encoding='win1251'");
    pqxx::result response = worker.exec("SELECT * FROM " + table_name);

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
    std::string table_name_drop;
    table_name_drop = '"' + table_name + '"';
    pqxx::result response = worker.exec("DROP TABLE IF EXISTS " + table_name);
    worker.commit();

    pqxx::work worker1(connectionObject);
    response = worker1.exec("CREATE TABLE " + table_name
        + " (id integer, x double precision, y double precision, z double precision, rad double precision, border integer, filled integer, emptied integer, path_to_border integer, way_to_border integer, neighbors_num integer[], neighbors_area double precision[]);");
    worker1.commit();
    connectionObject.close();
}
