/*
 * Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <iostream>

#include <spdlog/async.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace mc_rtc
{

namespace log
{

namespace details
{

inline spdlog::logger & success()
{
  static auto success = []() {
    if(spdlog::get("success"))
    {
      return spdlog::get("success");
    }
    auto success = spdlog::create_async_nb<spdlog::sinks::stdout_color_sink_mt>("success");
    success->set_pattern("%^[success]%$ %v");
    auto sink = static_cast<spdlog::sinks::stdout_color_sink_mt *>(success->sinks().back().get());
#ifndef WIN32
    sink->set_color(spdlog::level::info, "\033[01;32m"); // bold green
#else
    sink->set_color(spdlog::level::info, sink->BOLD | sink->GREEN);
#endif
    return success;
  }();
  return *success;
}

inline spdlog::logger & info()
{
  static auto info = []() {
    if(spdlog::get("info"))
    {
      return spdlog::get("info");
    }
    auto info = spdlog::create_async_nb<spdlog::sinks::stdout_color_sink_mt>("info");
    info->set_pattern("%^[info]%$ %v");
    auto sink = static_cast<spdlog::sinks::stdout_color_sink_mt *>(info->sinks().back().get());
#ifndef WIN32
    sink->set_color(spdlog::level::info, "\033[01;34m"); // bold cyan
#else
    sink->set_color(spdlog::level::info, sink->BOLD | sink->CYAN);
#endif
    return info;
  }();
  return *info;
}

inline spdlog::logger & cerr()
{
  static auto cerr = []() {
    if(spdlog::get("cerr"))
    {
      return spdlog::get("cerr");
    }
    auto cerr = spdlog::create_async_nb<spdlog::sinks::stderr_color_sink_mt>("cerr");
    cerr->set_pattern("[%^%l%$] %v");
    return cerr;
  }();
  return *cerr;
}

} // namespace details

template<typename ExceptionT, typename... Args>
void error_and_throw[[noreturn]](Args &&... args)
{
  auto message = fmt::format(std::forward<Args>(args)...);
  details::cerr().critical(message);
  throw ExceptionT(message);
}

template<typename... Args>
void error(Args &&... args)
{
  details::cerr().error(std::forward<Args>(args)...);
}

template<typename... Args>
void warning(Args &&... args)
{
  details::cerr().warn(std::forward<Args>(args)...);
}

template<typename... Args>
void info(Args &&... args)
{
  details::info().warn(std::forward<Args>(args)...);
}

template<typename... Args>
void success(Args &&... args)
{
  details::success().warn(std::forward<Args>(args)...);
}

} // namespace log

} // namespace mc_rtc

/*
#ifndef WIN32

 namespace mc_rtc
{

 constexpr auto OUT_NONE = "\033[00m";
 constexpr auto OUT_BLUE = "\033[01;34m";
 constexpr auto OUT_GREEN = "\033[01;32m";
 constexpr auto OUT_PURPLE = "\033[01;35m";
 constexpr auto OUT_RED = "\033[01;31m";

} // namespace mc_rtc

#  define LOG_ERROR(args) std::cerr << mc_rtc::OUT_RED << args << mc_rtc::OUT_NONE << std::endl;
#  define LOG_WARNING(args) std::cerr << mc_rtc::OUT_PURPLE << args << mc_rtc::OUT_NONE << std::endl;
#  define LOG_INFO(args) std::cout << mc_rtc::OUT_BLUE << args << mc_rtc::OUT_NONE << std::endl;
#  define LOG_SUCCESS(args) std::cout << mc_rtc::OUT_GREEN << args << mc_rtc::OUT_NONE << std::endl;

#else

#  include <windows.h>
 namespace mc_rtc
{
 static const HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
 constexpr auto OUT_NONE = 15;
 constexpr auto OUT_BLUE = 11;
 constexpr auto OUT_GREEN = 10;
 constexpr auto OUT_PURPLE = 13;
 constexpr auto OUT_RED = 12;
} // namespace mc_rtc

#  define LOG_ERROR(args)                                       \
    SetConsoleTextAttribute(mc_rtc::hConsole, mc_rtc::OUT_RED); \
    std::cerr << args << std::endl;                             \
    SetConsoleTextAttribute(mc_rtc::hConsole, mc_rtc::OUT_NONE);

#  define LOG_WARNING(args)                                        \
    SetConsoleTextAttribute(mc_rtc::hConsole, mc_rtc::OUT_PURPLE); \
    std::cerr << args << std::endl;                                \
    SetConsoleTextAttribute(mc_rtc::hConsole, mc_rtc::OUT_NONE);

#  define LOG_INFO(args)                                         \
    SetConsoleTextAttribute(mc_rtc::hConsole, mc_rtc::OUT_BLUE); \
    std::cout << args << std::endl;                              \
    SetConsoleTextAttribute(mc_rtc::hConsole, mc_rtc::OUT_NONE);

#  define LOG_SUCCESS(args)                                       \
    SetConsoleTextAttribute(mc_rtc::hConsole, mc_rtc::OUT_GREEN); \
    std::cout << args << std::endl;                               \
    SetConsoleTextAttribute(mc_rtc::hConsole, mc_rtc::OUT_NONE);

#endif

#define LOG_ERROR_AND_THROW(exception_type, args) \
  {                                               \
    std::stringstream strstrm;                    \
    strstrm << args;                              \
    LOG_ERROR(strstrm.str())                      \
    throw exception_type(strstrm.str());          \
  }
*/
