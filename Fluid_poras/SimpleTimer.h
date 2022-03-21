#pragma once
class SimpleTimer {
public:
	SimpleTimer() {
		start = std::chrono::high_resolution_clock::now();
	}
	~SimpleTimer() {
		end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> duration = end - start;
		std::cout << "Duratatioan: " << duration.count() << " s" << std::endl;
	}

private:
	std::chrono::time_point<std::chrono::steady_clock> start, end;
};
