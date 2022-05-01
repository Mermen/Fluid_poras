#pragma once
#include <utility>

#pragma once

class pore {
	double x_;
	double y_;
	double z_;
	double rad_;
	int border_;
	int filled_;
	int emptied_;
	int path_to_border_;
	int way_to_board_;
	std::vector<int> neighbors_num_;
	std::vector<double> neighbors_area_;
public: pore(const double x, const double y, const double z, const double rad, const int border, const int filled, const int emptied, const int path_to_border, const int way_to_board) {
	this->x_ = x;
	this->y_ = y;
	this->z_ = z;
	this->rad_ = rad;
	this->border_ = border;
	this->filled_ = filled;
	this->emptied_ = emptied;
	this->path_to_border_ = path_to_border;
	this->way_to_board_ = way_to_board;
}

	  double get_x() const
	  {
		  return this->x_;
	  }
	  double get_y() const
	  {
		  return this->y_;
	  }
	  double get_z() const
	  {
		  return this->z_;
	  }
	  double get_rad() const
	  {
		  return this->rad_;
	  }
	  int get_border() const
	  {
		  return this->border_;
	  }
	  int get_filled() const
	  {
		  return this->filled_;
	  }
	  int get_emptied() const
	  {
		  return this->emptied_;
	  }
	  int get_path_to_border() const
	  {
		  return this->path_to_border_;
	  }
	  int get_way_to_board() const
	  {
		  return this->way_to_board_;
	  }
	  int get_neighbor_num(const size_t i) const
	  {
		  return this->neighbors_num_[i];
	  }
	  std::vector<int> get_v_neighbor_num() {
		  return this->neighbors_num_;
	  }
	  double get_neighbor_area(const size_t i) const
	  {
		  return this->neighbors_area_[i];
	  }
	  std::vector<double> get_v_neighbor_area() {
		  return this->neighbors_area_;
	  }

	  void set_rad(const double rad) {
		  this->rad_ = rad;
	  }
	  void set_border(const int border) {
		  this->border_ = border;
	  }
	  void set_filled(const int filled) {
		  this->filled_ = filled;
	  }
	  void set_emptied(const int emptied) {
		  this->emptied_ = emptied;
	  }
	  void set_neighbor_num(const int k) {
		  this->neighbors_num_.push_back(k);
	  }
	  void set_neighbor_area(const double area) {
		  this->neighbors_area_.push_back(area);
	  }
	  void set_path_to_border(const int path_to_border) {
		  this->path_to_border_ = path_to_border;
	  }
	  void set_way_to_board(const int way_to_board) {
		  this->way_to_board_ = way_to_board;
	  }
	  void set_neighbors_num(std::vector<int> neighbors_num) {
		  this->neighbors_num_ = std::move(neighbors_num);
	  }
	  void set_neighbors_area(std::vector<double> neighbors_area) {
		  this->neighbors_area_ = std::move(neighbors_area);
	  }

	  int find_k(const int k) const
	  {
		  for (const int i : this->neighbors_num_)
		  {
			  if (i == k) {
				  return 1;

			  }
		  }
		  return 0;

	  }

	  int road_to_board(std::vector<pore>& pores) {

		  if (this->border_ || this->way_to_board_)
			  return 1;
		  else if (this->emptied_|| this->path_to_border_)
			  return 0;
		  else {
			  this->path_to_border_ = 1;
			  int road = 0;
			  for (size_t i = 0; i < this->neighbors_num_.size(); i++) {
				  road += pores[neighbors_num_[i]].road_to_board(pores);
			  }
			  if (road) {
				  this->way_to_board_ = 1;
				  return 1;
			  }
			  else {
				  return 0;
			  }
		  }
	  }

	  int filled_now(const std::vector<pore>& pores) const
	  {
		  int ret = 0;
		  for (const int i : this->neighbors_num_)
		  {
			  if (pores[i].get_filled()) {
				  ret++;
			  }
		  }
		  return ret;
	  }
};
