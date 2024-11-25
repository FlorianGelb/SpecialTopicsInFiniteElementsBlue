//
// Created by Florian on 25.11.2024.
//

#ifndef LAPLACIAN_TIMEHANDLER_H
#define LAPLACIAN_TIMEHANDLER_H


class TimeHandler
{
public:
  TimeHandler(const double time_end, const double delta_t)
    : timestep(0)
    , time_current(0.0)
    , time_end(time_end)
    , delta_t(delta_t)
  {}

  virtual ~TimeHandler() = default;

  double current() const
  {
    return time_current;
  }
  double end() const
  {
    return time_end;
  }
  double get_delta_t() const
  {
    return delta_t;
  }
  unsigned int get_timestep() const
  {
    return timestep;
  }
  void increment()
  {
    time_current += delta_t;
    ++timestep;
  }

private:
  unsigned int timestep;
  double       time_current;
  const double time_end;
  const double delta_t;
};

#endif // LAPLACIAN_TIMEHANDLER_H
