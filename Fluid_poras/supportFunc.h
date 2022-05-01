#pragma once
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

inline double linear(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& k, const std::vector<double>& b, const double X, const size_t n) {
	// за пределами точек
	if (X <= x[0])
		return y[0];
	else if (X >= x[n - 1])
		return y[n - 1];

	int left = 0;
	int right = n;

	while (true) {
		if (const int middle = (left + right) / 2; X < x[middle]) {
			right = middle;
		}
		else if (X > x[middle]) {
			left = middle;
		}
		if ((left + 1) == right) {
			return k[left] * X + b[left];
		}

	}
}

inline std::string to_string_pore(pore pore) {
	std::string str;
	std::string str_num;
	std::string str_area;
	std::ostringstream o_str_num;
	std::ostringstream o_str_area;
	std::vector<int> vec_num;
	std::vector<double> vec_area;
	vec_num = pore.get_v_neighbor_num();
	vec_area = pore.get_v_neighbor_area();

	str_num += '\'';
	str_num += '{';
	if (!vec_num.empty()) {
		// Convert all but the last element to avoid a trailing ","
		std::copy(vec_num.begin(), vec_num.end() - 1,
			std::ostream_iterator<int>(o_str_num, ","));

		// Now add the last element with no delimiter
		o_str_num << vec_num.back();
	}
	str_num += o_str_num.str();
	str_num += '}';
	str_num += '\'';

	str_area += '\'';
	str_area += '{';
	if (!vec_area.empty()) {
		// Convert all but the last element to avoid a trailing ","
		std::copy(vec_area.begin(), vec_area.end() - 1,
			std::ostream_iterator<double>(o_str_area, ","));

		// Now add the last element with no delimiter
		o_str_area << vec_area.back();
	}
	str_area += o_str_area.str();
	str_area += '}';
	str_area += '\'';

	str = std::to_string(pore.get_x()) + "," + std::to_string(pore.get_y()) + "," + std::to_string(pore.get_z()) + "," + std::to_string(pore.get_rad()) + ","
		+ std::to_string(pore.get_border()) + "," + std::to_string(pore.get_filled()) + "," + std::to_string(pore.get_emptied()) + ","
		+ std::to_string(pore.get_path_to_border()) + "," + std::to_string(pore.get_way_to_board()) + ","
		+ str_num + "," + str_area;
	return str;
}

template <typename T>
bool find_bool(std::vector<T>& vec, T elem) {
	auto result = std::find(begin(vec), end(vec), elem);
	return (result != std::end(vec));
}