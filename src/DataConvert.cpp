#include "DataConvert.h"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

int64_t DataConvert::safeStoll(const std::string& str)
{
	try {
        if (str.empty()) {
            return 0.0; // ����������Ҫ��Ĭ��ֵ������ NaN
        }
        std::string processedStr = str;
        if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
            processedStr = processedStr.substr(1, processedStr.size() - 2);
        }
        std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
            [](unsigned char c) { return std::tolower(c); });
        std::istringstream iss(processedStr);
        if (str == "null") {
            return std::numeric_limits<int64_t>::min(); // ʹ�� INT_MIN ��ʾ������
        }

        long long num;
        iss >> num;
        // ����Ƿ񳬳� double �ı�ʾ��Χ
        if (num > std::numeric_limits<int64_t>::max() || num < -std::numeric_limits<int64_t>::min()) {
            std::cerr << "Out of range: Value exceeds double representation for string: " << str << std::endl;
            return std::numeric_limits<int64_t>::min(); // ���� INT_MIN
        }
        int64_t value = static_cast<int64_t>(num);
        return value;
	}
	catch (const std::invalid_argument& e) {
		// ������Ч�����쳣 (���磬�ַ���������Ч������)
		std::cerr << "Invalid argument: " << e.what() << " for string: " << str << std::endl;
		return  std::numeric_limits<int64_t>::min(); // ����������Ҫ��Ĭ��ֵ
	}
	catch (const std::out_of_range& e) {
		// ��������Χ�쳣 (���磬����̫���̫С)
		std::cerr << "Out of range: " << e.what() << " for string: " << str << std::endl;
		return  std::numeric_limits<int64_t>::min(); // ����������Ҫ��Ĭ��ֵ
	}
}

int16_t DataConvert::safeStoi16_t(const std::string& str)
{
    try {
        if (str.empty()) {
            return 0.0; // ����������Ҫ��Ĭ��ֵ������ NaN
        }
        std::string processedStr = str;
        if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
            processedStr = processedStr.substr(1, processedStr.size() - 2);
        }
        std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
            [](unsigned char c) { return std::tolower(c); });
        std::istringstream iss(processedStr);
        if (str == "null") {
            return std::numeric_limits<int16_t>::min(); // ʹ�� INT_MIN ��ʾ������
        }

        long long num;
        iss >> num;
        // ����Ƿ񳬳� double �ı�ʾ��Χ
        if (num > std::numeric_limits<int16_t>::max() || num < std::numeric_limits<int16_t>::min()) {
            std::cerr << "Out of range: Value exceeds double representation for string: " << str << std::endl;
            return std::numeric_limits<int16_t>::min(); // ���� NaN
        }
        int16_t value = static_cast<int16_t>(num);
        return value;
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << " for string: " << str << std::endl;
        return std::numeric_limits<int16_t>::min(); // ���� NaN
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Out of range: " << e.what() << " for string: " << str << std::endl;
        return std::numeric_limits<int16_t>::min(); // ���� NaN
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
        return std::numeric_limits<int8_t>::min(); // ʹ�� INT_MIN ��ʾ������
    }

    long long num; // ʹ��long long�����м������
    iss >> num;

    if (num < std::numeric_limits<int8_t>::min() || num > std::numeric_limits<int8_t>::max()) {
        // �����������������ַ����������֣����߳���int8_t�ķ�Χ
        // �����׳��쳣������Ĭ��ֵ����������������ʽ
        // ����ѡ���׳��쳣
        throw std::runtime_error("Invalid string or out of range for int8_t");
    }

    return static_cast<int8_t>(num);
}

double DataConvert::safeStod(const std::string& str) {
    try {
        if (str.empty()) {
            return 0.0; // ����������Ҫ��Ĭ��ֵ������ NaN
        }

        std::string processedStr = str;
        if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
            processedStr = processedStr.substr(1, processedStr.size() - 2);
        }
        std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
            [](unsigned char c) { return std::tolower(c); });
        std::istringstream iss(processedStr);
        if (str == "null") {
            return -std::numeric_limits<double>::infinity(); // ʹ�� nan ��ʾ������
        }

        double num; // ʹ��long long�����м������
        iss >> num;
        if (num < -std::numeric_limits<double>::infinity() || num > std::numeric_limits<double>::infinity()) {
            throw std::runtime_error("Invalid string or out of range for double");
        }

        return static_cast<double>(num);
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << " for string: " << str << std::endl;
        return  -std::numeric_limits<float>::infinity(); // ���� NaN
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Out of range: " << e.what() << " for string: " << str << std::endl;
        return  -std::numeric_limits<float>::infinity(); // ���� NaN
    }
}

float DataConvert::safeStof(const std::string& str) {
    try {
        if (str.empty()) {
            return std::numeric_limits<float>::quiet_NaN();  // ���� NaN ������Ĭ��ֵ
        }

        std::string processedStr = str;
        if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
            processedStr = processedStr.substr(1, processedStr.size() - 2);
        }
        std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
            [](unsigned char c) { return std::tolower(c); });
        std::istringstream iss(processedStr);
        if (str == "null") {
            return -std::numeric_limits<float>::infinity(); // ʹ�� nan ��ʾ������
        }
        float num; // ʹ��long long�����м������
        iss >> num;
        if (num < -std::numeric_limits<float>::infinity() || num > std::numeric_limits<float>::infinity()) {
            // �����������������ַ����������֣����߳���int8_t�ķ�Χ
            // �����׳��쳣������Ĭ��ֵ����������������ʽ
            // ����ѡ���׳��쳣
            throw std::runtime_error("Invalid string or out of range for float");
        }
        return static_cast<float>(num);
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << " for string: " << str << std::endl;
        return  -std::numeric_limits<float>::infinity(); // ���� NaN
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Out of range: " << e.what() << " for string: " << str << std::endl;
        return  -std::numeric_limits<float>::infinity(); // ���� NaN
    }
}

bool DataConvert::safeStob(const std::string& str)
{
    std::string processedStr = str;

    // ��鲢ȥ����β��˫����
    if (processedStr.size() >= 2 && processedStr.front() == '"' && processedStr.back() == '"') {
        processedStr = processedStr.substr(1, processedStr.size() - 2);
    }

    // ת��ΪСд
    std::transform(processedStr.begin(), processedStr.end(), processedStr.begin(),
        [](unsigned char c) { return std::tolower(c); });

    return processedStr == "true";
}
