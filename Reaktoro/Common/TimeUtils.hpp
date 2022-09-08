
#pragma once

// C++ includes
#include <chrono>

namespace Reaktoro {

using Time = std::chrono::time_point<std::chrono::high_resolution_clock>;

using Duration = std::chrono::duration<double>;

/// Return the time point now
/// @see elapsed
auto time() -> Time;

/// Return the elapsed time between two time points (in units of s)
/// @param end The end time point
/// @param end The begin time point
/// @return The elapsed time between *end* and *begin* in seconds
auto elapsed(const Time& end, const Time& begin) -> double;

/// Return the elapsed time between a time point and now (in units of s)
/// @param end The begin time point
/// @return The elapsed time between now and *begin* in seconds
auto elapsed(const Time& begin) -> double;

/// Used for measuring elapsed time since object creation.
class Stopwatch
{
public:
    /// Construct a Stopwatch object and start measuring time.
    /// @param reading The variable to store the final measured elapsed time
    Stopwatch(double& reading) : reading(reading), tstart(time()) {}

    /// Destroy this Stopwatch object and stop measuring time.
    ~Stopwatch() { stop(); }

    /// Start measuring time.
    auto start() -> void { tstart = time(); }

    /// Stop measuring time.
    auto stop() -> void { reading = elapsed(tstart); }

private:
    /// The auxiliary time variable marking the start of timing.
    Time tstart;

    /// The reference to the variable that will store the final measured elapsed time.
    double& reading;
};

} // namespace Reaktoro
