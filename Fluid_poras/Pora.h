#pragma once

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

	  void set_Rad(double Rad) {
		  this->Rad = Rad;
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
		  for (int i = 0; i < this->neighbors_num.size(); i++) {
			  if (this->neighbors_num[i] == k) {
				  return 1;

			  }
		  }
		  return 0;

	  }

	  int road_to_board(std::vector<Pora>& poras) {

		  if (this->border) {
			  return 1;
		  }
		  else if (this->emptied) {
			  return 0;
		  }
		  else if (this->way_to_board) {
			  return 1;
		  }
		  else if (this->path_to_border) {
			  return 0;
		  }
		  else {
			  this->path_to_border = 1;
			  int road = 0;
			  for (int i = 0; i < this->neighbors_num.size(); i++) {
				  road += poras[neighbors_num[i]].road_to_board(poras);
			  }
			  if (road) {
				  this->way_to_board = 1;
				  return 1;
			  }
			  else {
				  return 0;
			  }
		  }
	  }

	  int filled_now(std::vector<Pora>& poras) {
		  int ret = 0;
		  for (size_t i = 0; i < this->neighbors_num.size(); i++) {
			  if (poras[this->neighbors_num[i]].get_filled()) {
				  ret++;
			  }
		  }
		  return ret;
	  }
};
