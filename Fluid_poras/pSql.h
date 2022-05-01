#pragma once
#include <utility>

#pragma once
inline void insert_into_table(const std::vector<pore>& pores, const std::string& table_name, const std::string& connection_string) {
	pqxx::connection connection_object(connection_string.c_str());
	pqxx::work worker(connection_object);
	worker.exec("set client_encoding='win1251'");
	for (size_t i = 0; i < pores.size(); i++) {
		worker.exec("INSERT INTO " + table_name + " values(" + std::to_string(i) + ","
			+ to_string_pore(pores[i]) + ");");
	}
	worker.commit();
	connection_object.close();
}

inline void insert_into_table_graph(const std::string& graph_id, std::vector<double> volume_filled_graph, std::vector<double> pressure_graph,
                                    const std::string& table_name, const std::string& connection_string) {
	std::string str_volume;
	std::string str_pressure;
	std::ostringstream o_str_volume;
	std::ostringstream o_str_pressure;
	std::vector<double> vec_volume;
	std::vector<double> vec_pressure;
	vec_volume = std::move(volume_filled_graph);
	vec_pressure = std::move(pressure_graph);

	str_volume += '\'';
	str_volume += '{';
	if (!vec_volume.empty()) {
		// Convert all but the last element to avoid a trailing ","
		std::copy(vec_volume.begin(), vec_volume.end() - 1,
			std::ostream_iterator<double>(o_str_volume, ","));

		// Now add the last element with no delimiter
		o_str_volume << vec_volume.back();
	}
	str_volume += o_str_volume.str();
	str_volume += '}';
	str_volume += '\'';

	str_pressure += '\'';
	str_pressure += '{';
	if (!vec_pressure.empty()) {
		// Convert all but the last element to avoid a trailing ","
		std::copy(vec_pressure.begin(), vec_pressure.end() - 1,
			std::ostream_iterator<double>(o_str_pressure, ","));

		// Now add the last element with no delimiter
		o_str_pressure << vec_pressure.back();
	}
	str_pressure += o_str_pressure.str();
	str_pressure += '}';
	str_pressure += '\'';
	/*
	std::cout << table_name << std::endl;
	std::cout << graph_id << std::endl;
	std::cout << str_volume << std::endl;
	std::cout << str_pressure << std::endl;
	*/

	pqxx::connection connection_object(connection_string.c_str());
	pqxx::work worker(connection_object);
	worker.exec("set client_encoding='win1251'");
	worker.exec("DELETE FROM " + table_name + " WHERE id=" + graph_id + ";");
	worker.exec("INSERT INTO " + table_name + " values(" + graph_id + ","
		+ str_volume + "," + str_pressure + ",'url');");
	worker.commit();
	connection_object.close();
}

inline void insert_into_table_graph_li(const std::string& graph_id, std::vector<int> porous_neighbours, const std::string&
                                       connection_string) {
	std::string str_neighbours;
	std::ostringstream o_str_neighbours;
	std::vector<int> vec_neighbours = std::move(porous_neighbours);

	str_neighbours += '\'';
	str_neighbours += '{';
	if (!vec_neighbours.empty()) {
		// Convert all but the last element to avoid a trailing ","
		std::copy(vec_neighbours.begin(), vec_neighbours.end() - 1,
			std::ostream_iterator<double>(o_str_neighbours, ","));

		// Now add the last element with no delimiter
		o_str_neighbours << vec_neighbours.back();
	}
	str_neighbours += o_str_neighbours.str();
	str_neighbours += '}';
	str_neighbours += '\'';


	pqxx::connection connectionObject(connection_string.c_str());
	pqxx::work worker(connectionObject);
	worker.exec("set client_encoding='win1251'");
	worker.exec("DELETE FROM graph_li WHERE id=" + graph_id + ";");
	worker.exec("INSERT INTO graph_li values(" + graph_id + ","
		+ str_neighbours + ");");
	worker.commit();
	connectionObject.close();
}

inline void select_from_table(std::vector<pore>& pores, const std::string& table_name, const std::string& connection_string) {
	std::string str;
	std::string tmp;
	pqxx::connection connection_object(connection_string.c_str());
	pqxx::work worker(connection_object);
	worker.exec("set client_encoding='win1251'");
	pqxx::result response = worker.exec("SELECT * FROM " + table_name+" ORDER BY id");

	for (int i = 0; i < response.size(); i++) {
		pores.emplace_back(std::stod(to_string(response[i][1])), std::stod(to_string(response[i][2])), std::stod(to_string(response[i][3])), std::stod(to_string(response[i][4])),
		                   std::stoi(to_string(response[i][5])), std::stoi(to_string(response[i][6])), std::stoi(to_string(response[i][7])), std::stoi(to_string(response[i][8])), std::stoi(to_string(response[i][9])));
		std::stringstream ss0(to_string(response[i][10]).substr(1, to_string(response[i][10]).size() - 2));
		std::vector<int> neighbors_num;
		while (std::getline(ss0, tmp, ',')) {
			neighbors_num.push_back(std::stoi(tmp));
		}
		pores[i].set_neighbors_num(neighbors_num);

		std::stringstream ss1(to_string(response[i][11]).substr(1, to_string(response[i][11]).size() - 2));
		std::vector<double> neighbors_area;
		while (std::getline(ss1, tmp, ',')) {
			neighbors_area.push_back(std::stod(tmp));
		}
		pores[i].set_neighbors_area(neighbors_area);

		/*
		str += to_string(response[i][0]) + " ";
		for (int j = 1; j < 5; j++)
		{
			str += to_string(response[i][j]) + " ";
		}
		for (int j = 5; j < 10; j++)
		{
			str += to_string(response[i][j]) + " ";
		}

		std::stringstream ss0(to_string(response[i][10]).substr(1, to_string(response[i][10]).size() - 2));
		std::vector<int> neighbors_num;
		while (std::getline(ss0, tmp, ',')) {
			neighbors_num.push_back(std::stoi(tmp));
		}
		for (int j = 0; j < neighbors_num.size(); j++)
		{
			str += std::to_string(neighbors_num[j])+" ";
		}

		std::stringstream ss1(to_string(response[i][11]).substr(1, to_string(response[i][11]).size() - 2));
		std::vector<double> neighbors_area;
		while (std::getline(ss1, tmp, ',')) {
			neighbors_area.push_back(std::stod(tmp));
		}
		for (int j = 0; j < neighbors_area.size(); j++)
		{
			str += std::to_string(neighbors_area[j]) + " ";
		}
		std::cout << str << std::endl;
		str = "";
		*/
	}
	worker.commit();
	connection_object.close();
}

inline void create_table(const std::string& table_name, const std::string& connection_string) {
	std::string str;
	std::string tmp;
	pqxx::connection connection_object(connection_string.c_str());
	pqxx::work worker(connection_object);
	worker.exec("set client_encoding='win1251'");
	worker.exec("DROP TABLE IF EXISTS " + table_name);
	worker.commit();

	pqxx::work worker1(connection_object);
	worker1.exec("CREATE TABLE " + table_name
		+ " (id integer, x double precision, y double precision, z double precision, rad double precision, border integer, filled integer, emptied integer, path_to_border integer, way_to_border integer, neighbors_num integer[], neighbors_area double precision[]);");
	worker1.commit();
	connection_object.close();
}
