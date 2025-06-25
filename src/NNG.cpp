#include "NNG.h"

NNG::NNG() :map()
{

}

NNG::~NNG()
{

}

void NNG::push(Edge* e)
{
	this->map.insert(std::make_pair(e->weight, e));
}

Edge* NNG::pop()
{
	if (!this->map.empty())
	{
		auto it = this->map.begin();
		Edge* min_edge = it->second;
		this->map.erase(it);
		return min_edge;
	}
	else
	{
		return nullptr;
	}

}

bool NNG::find(Edge* ptrToFind)
{
	auto range = this->map.equal_range(ptrToFind->weight);
	for (auto it = range.first; it != range.second; ++it) {
		if (it->second == ptrToFind) {
			/*std::cout << "A* found in the unordered multimap." << std::endl;*/
			return true;
		}
	}
	return false;
}

void NNG::remove(Edge* edge_to_remove)
{
	for (auto it = map.begin(); it != map.end();) {
		if (it->second == edge_to_remove) 
		{
			this->map.erase(it); // edge的地址都是唯一的，直接删除元素
			break;
		}
		else 
		{
			++it;
		}
	}
}

int NNG::size()
{
	return this->map.size();
}

bool NNG::empty()
{
	return this->map.empty();
} 

void NNG::clear()
{
	this->map.clear();
}

void NNG::show_map()
{
	std::cout << "NNG: size="<<this->size()<<std::endl;
	for (const auto& entry : this->map) {
		std::cout <<"(ID=" << entry.second->ID << ", weight=" << entry.first << ") ";
	}
	std::cout << std::endl;
}