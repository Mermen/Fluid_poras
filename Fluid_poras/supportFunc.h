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

std::string to_string_pora(Pora pora) {
    std::string str;
    std::string str_num = "";
    std::string str_area = "";
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
            std::ostream_iterator<double>(ostr_area, ","));

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

bool find_bool(std::vector<int>& vec, int elem) {
    auto result = std::find(begin(vec), end(vec), elem);
    if (result != std::end(vec)) {
        return 1;
    }
    else {
        return 0;
    }
}