#pragma once
class simple_timer {
public:
	simple_timer() {
		start_ = std::chrono::high_resolution_clock::now();
	}
	~simple_timer() {
		end_ = std::chrono::high_resolution_clock::now();
		const std::chrono::duration<float> duration = end_ - start_;
		std::cout << "Duration: " << duration.count() << " s" << std::endl;
	}

private:
	std::chrono::time_point<std::chrono::steady_clock> start_, end_;
};
