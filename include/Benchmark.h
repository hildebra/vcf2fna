#pragma once

#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>

class Benchmark {
public:
	explicit Benchmark(std::string label) :
		name(std::move(label)), timeSumMicroseconds(0),
		startTime(Clock::now()), lastTime(startTime) {}

	std::uint64_t now_total_secs() const {
		const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(Clock::now() - startTime);
		return elapsed.count() < 0 ? 0U : static_cast<std::uint64_t>(elapsed.count());
	}

	float now_interval_secs(int decimals = 0) {
		const auto now = Clock::now();
		const std::chrono::duration<double> elapsed = now - lastTime;
		lastTime = now;
		double seconds = elapsed.count();
		if (decimals > 0) {
			const double multiplier = std::pow(10.0, static_cast<double>(decimals));
			seconds = std::ceil(seconds * multiplier) / multiplier;
		}
		return static_cast<float>(seconds);
	}

	void start() {
		startTime = Clock::now();
		lastTime = startTime;
	}

	void stop() {
		const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - startTime);
		timeSumMicroseconds += elapsed.count();
	}

	void now_total_time(std::ostream& output = std::cout) {
		const auto now = Clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(now - lastTime);
		lastTime = now;
		printDuration(elapsed.count(), output);
	}

	void printResults(std::ostream& output = std::cout) const {
		printDuration(timeSumMicroseconds, output);
	}

private:
	using Clock = std::chrono::steady_clock;

	void printDuration(std::int64_t microseconds, std::ostream& output) const {
		if (microseconds < 0) microseconds = 0;
		const std::int64_t hours = microseconds / 3'600'000'000LL;
		microseconds %= 3'600'000'000LL;
		const std::int64_t minutes = microseconds / 60'000'000LL;
		microseconds %= 60'000'000LL;
		const std::int64_t seconds = microseconds / 1'000'000LL;
		microseconds %= 1'000'000LL;
		const std::int64_t milliseconds = microseconds / 1'000LL;

		output << name << ": ";
		if (hours != 0) output << hours << "h ";
		if (minutes != 0) output << minutes << "m ";
		if (seconds != 0) output << seconds << "s ";
		if (milliseconds != 0 || (hours == 0 && minutes == 0 && seconds == 0)) {
			output << milliseconds << "ms";
		}
		output << '\n';
	}

	std::string name;
	std::int64_t timeSumMicroseconds;
	Clock::time_point startTime;
	Clock::time_point lastTime;
};
