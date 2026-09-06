#pragma once
#include "Data_structs.h"
#include <vector>
#include <array>

std::array<double, 144> get_element_K(Material mat, double dx, double dy);
void check_capacities(const std::vector<Element>& els, const std::vector<double>& U, Material mat, double dx, double dy, double& mmx, double& mmy);