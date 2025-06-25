#include "DataConvert.h"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

int64_t DataConvert::safeStoll(const std::string& str)
{
	try {
        if (str.empty()) {
            return 0.0; // 或其他你需要的默认值，例如 NaN
        }
        std::string processedStr = str;
        if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
            processedStr = processedStr.substr(1, processedStr.size() - 2);
        }
        std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
            [](unsigned char c) { return std::tolower(c); });
        std::istringstream iss(processedStr);
        if (str == "null") {
            return std::numeric_limits<int64_t>::min(); // 使用 INT_MIN 表示空数据
        }

        long long num;
        iss >> num;
        // 检查是否超出 double 的表示范围
        if (num > std::numeric_limits<int64_t>::max() || num < -std::numeric_limits<int64_t>::min()) {
            std::cerr << "Out of range: Value exceeds double representation for string: " << str << std::endl;
            return std::numeric_limits<int64_t>::min(); // 返回 INT_MIN
        }
        int64_t value = static_cast<int64_t>(num);
        return value;
	}
	catch (const std::invalid_argument& e) {
		// 处理无效参数异常 (例如，字符串不是有效的数字)
		std::cerr << "Invalid argument: " << e.what() << " for string: " << str << std::endl;
		return  std::numeric_limits<int64_t>::min(); // 或其他你需要的默认值
	}
	catch (const std::out_of_range& e) {
		// 处理超出范围异常 (例如，数字太大或太小)
		std::cerr << "Out of range: " << e.what() << " for string: " << str << std::endl;
		return  std::numeric_limits<int64_t>::min(); // 或其他你需要的默认值
	}
}

int16_t DataConvert::safeStoi16_t(const std::string& str)
{
    try {
        if (str.empty()) {
            return 0.0; // 或其他你需要的默认值，例如 NaN
        }
        std::string processedStr = str;
        if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
            processedStr = processedStr.substr(1, processedStr.size() - 2);
        }
        std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
            [](unsigned char c) { return std::tolower(c); });
        std::istringstream iss(processedStr);
        if (str == "null") {
            return std::numeric_limits<int16_t>::min(); // 使用 INT_MIN 表示空数据
        }

        long long num;
        iss >> num;
        // 检查是否超出 double 的表示范围
        if (num > std::numeric_limits<int16_t>::max() || num < std::numeric_limits<int16_t>::min()) {
            std::cerr << "Out of range: Value exceeds double representation for string: " << str << std::endl;
            return std::numeric_limits<int16_t>::min(); // 返回 NaN
        }
        int16_t value = static_cast<int16_t>(num);
        return value;
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << " for string: " << str << std::endl;
        return std::numeric_limits<int16_t>::min(); // 返回 NaN
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Out of range: " << e.what() << " for string: " << str << std::endl;
        return std::numeric_limits<int16_t>::min(); // 返回 NaN
    }
}

int8_t DataConvert::safeStoi8_t(const std::string& str)
{
    std::string processedStr = str;
    if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
        processedStr = processedStr.substr(1, processedStr.size() - 2);
    }
    std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
        [](unsigned char c) { return std::tolower(c); });
    std::istringstream iss(processedStr);
    if (str == "null") {
        return std::numeric_limits<int8_t>::min(); // 使用 INT_MIN 表示空数据
    }

    long long num; // 使用long long避免中间结果溢出
    iss >> num;

    if (num < std::numeric_limits<int8_t>::min() || num > std::numeric_limits<int8_t>::max()) {
        // 处理错误情况，例如字符串不是数字，或者超出int8_t的范围
        // 可以抛出异常，返回默认值，或者其他错误处理方式
        // 这里选择抛出异常
        throw std::runtime_error("Invalid string or out of range for int8_t");
    }

    return static_cast<int8_t>(num);
}

double DataConvert::safeStod(const std::string& str) {
    try {
        if (str.empty()) {
            return 0.0; // 或其他你需要的默认值，例如 NaN
        }

        std::string processedStr = str;
        if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
            processedStr = processedStr.substr(1, processedStr.size() - 2);
        }
        std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
            [](unsigned char c) { return std::tolower(c); });
        std::istringstream iss(processedStr);
        if (str == "null") {
            return -std::numeric_limits<double>::infinity(); // 使用 nan 表示空数据
        }

        double num; // 使用long long避免中间结果溢出
        iss >> num;
        if (num < -std::numeric_limits<double>::infinity() || num > std::numeric_limits<double>::infinity()) {
            throw std::runtime_error("Invalid string or out of range for double");
        }

        return static_cast<double>(num);
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << " for string: " << str << std::endl;
        return  -std::numeric_limits<float>::infinity(); // 返回 NaN
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Out of range: " << e.what() << " for string: " << str << std::endl;
        return  -std::numeric_limits<float>::infinity(); // 返回 NaN
    }
}

float DataConvert::safeStof(const std::string& str) {
    try {
        if (str.empty()) {
            return std::numeric_limits<float>::quiet_NaN();  // 返回 NaN 或其他默认值
        }

        std::string processedStr = str;
        if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
            processedStr = processedStr.substr(1, processedStr.size() - 2);
        }
        std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
            [](unsigned char c) { return std::tolower(c); });
        std::istringstream iss(processedStr);
        if (str == "null") {
            return -std::numeric_limits<float>::infinity(); // 使用 nan 表示空数据
        }
        float num; // 使用long long避免中间结果溢出
        iss >> num;
        if (num < -std::numeric_limits<float>::infinity() || num > std::numeric_limits<float>::infinity()) {
            // 处理错误情况，例如字符串不是数字，或者超出int8_t的范围
            // 可以抛出异常，返回默认值，或者其他错误处理方式
            // 这里选择抛出异常
            throw std::runtime_error("Invalid string or out of range for float");
        }
        return static_cast<float>(num);
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << " for string: " << str << std::endl;
        return  -std::numeric_limits<float>::infinity(); // 返回 NaN
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Out of range: " << e.what() << " for string: " << str << std::endl;
        return  -std::numeric_limits<float>::infinity(); // 返回 NaN
    }
}

bool DataConvert::safeStob(const std::string& str)
{
    std::string processedStr = str;

    // 检查并去除首尾的双引号
    if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
        processedStr = processedStr.substr(1, processedStr.size() - 2);
    }

    // 转换为小写
    std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
        [](unsigned char c) { return std::tolower(c); });

    return processedStr == "true";
}
