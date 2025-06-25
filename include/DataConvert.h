#pragma once
#include <string>
class DataConvert
{
public:
	static int64_t safeStoll(const std::string& str);
	static int16_t safeStoi16_t(const std::string& str);
	static int8_t safeStoi8_t(const std::string& str);
	static double safeStod(const std::string& str);
	static float safeStof(const std::string& str);
	static bool safeStob(const std::string& str);
};

