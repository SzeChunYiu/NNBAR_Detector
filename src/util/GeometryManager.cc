// ============================================================================
// GeometryManager.cc
// Implementation of geometry lookup and management
// ============================================================================

#include "util/GeometryManager.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace nnbar {

// ============================================================================
// Singleton instance
// ============================================================================

GeometryManager& GeometryManager::Instance() {
    static GeometryManager instance;
    return instance;
}

// ============================================================================
// Initialization
// ============================================================================

void GeometryManager::Initialize() {
    if (m_initialized) {
        std::cout << "GeometryManager: Already initialized, clearing..." << std::endl;
        Clear();
    }

    std::cout << "GeometryManager: Initializing geometry lookup tables" << std::endl;
    m_initialized = true;
}

void GeometryManager::Clear() {
    m_scintModules.clear();
    m_scintBars.clear();
    m_leadGlassBlocks.clear();
    m_tpcModules.clear();
    m_beampipeSections.clear();
    m_scintBarIndex.clear();
    m_initialized = false;
}

// ============================================================================
// Registration methods
// ============================================================================

void GeometryManager::RegisterScintillatorModule(int32_t copyNumber,
                                                  double x, double y, double z,
                                                  int32_t surfaceIndex,
                                                  double rotationAngle) {
    VolumeInfo info;
    info.copyNumber = copyNumber;
    info.name = "ScintPV";
    info.type = "scint_module";
    info.x = x;
    info.y = y;
    info.z = z;
    info.surfaceIndex = surfaceIndex;
    info.rotZ = rotationAngle;
    info.energyDeposit = 0.0;
    info.nOpticalPhotons = 0;

    m_scintModules.push_back(info);
}

void GeometryManager::RegisterScintillatorBar(int32_t moduleCopyNumber,
                                               int32_t barCopyNumber,
                                               const std::string& barType,
                                               double localX, double localY, double localZ,
                                               double barSizeX, double barSizeY, double barSizeZ) {
    // Find the module
    const VolumeInfo* moduleInfo = nullptr;
    for (const auto& mod : m_scintModules) {
        if (mod.copyNumber == moduleCopyNumber) {
            moduleInfo = &mod;
            break;
        }
    }

    ScintillatorInfo bar;
    bar.copyNumber = moduleCopyNumber;
    bar.barCopyNumber = barCopyNumber;
    bar.barType = barType;
    bar.localX = localX;
    bar.localY = localY;
    bar.localZ = localZ;
    bar.barX = barSizeX;
    bar.barY = barSizeY;
    bar.barZ = barSizeZ;
    bar.energyDeposit = 0.0;
    bar.nPhotons = 0;

    if (moduleInfo) {
        bar.moduleX = moduleInfo->x;
        bar.moduleY = moduleInfo->y;
        bar.moduleZ = moduleInfo->z;
        bar.surfaceIndex = moduleInfo->surfaceIndex;
        bar.rotationAngle = moduleInfo->rotZ;

        // Calculate global position with rotation
        double rotRad = bar.rotationAngle * M_PI / 180.0;
        double cosR = std::cos(rotRad);
        double sinR = std::sin(rotRad);

        // Apply rotation around Z
        bar.globalX = bar.moduleX + (localX * cosR - localY * sinR);
        bar.globalY = bar.moduleY + (localX * sinR + localY * cosR);
        bar.globalZ = bar.moduleZ + localZ;
    }

    // Store and index
    size_t idx = m_scintBars.size();
    m_scintBars.push_back(bar);
    m_scintBarIndex[{moduleCopyNumber, barCopyNumber}] = idx;
}

void GeometryManager::RegisterLeadGlassBlock(int32_t copyNumber,
                                              double x, double y, double z,
                                              double rotX, double rotZ,
                                              int32_t surfaceIndex) {
    LeadGlassInfo info;
    info.copyNumber = copyNumber;
    info.csvX = x;
    info.csvY = y;
    info.csvZ = z;
    info.csvRotX = rotX;
    info.csvRotZ = rotZ;
    info.surfaceIndex = surfaceIndex;
    info.energyDeposit = 0.0;
    info.nPhotons = 0;
    info.nPMTHits = 0;

    // Calculate global position based on surface
    // Surface rotation: 0=top, 1=right(90), 2=bottom(180), 3=left(270)
    if (surfaceIndex < 4) {
        // Side surfaces
        double angle = surfaceIndex * 90.0;
        double angleRad = angle * M_PI / 180.0;
        double cosA = std::cos(angleRad);
        double sinA = std::sin(angleRad);

        // Original position is (x, y+offset, z) where offset accounts for radial position
        // Rotate around Z axis by surface angle
        info.globalX = x * 10.0 * cosA - (y * 10.0 + 15.0) * sinA;  // cm to mm, offset added
        info.globalY = x * 10.0 * sinA + (y * 10.0 + 15.0) * cosA;
        info.globalZ = z * 10.0;  // cm to mm
        info.rotationAngle = angle;
    } else {
        // Front (4) or Back (5) surfaces
        info.globalX = x * 10.0;
        info.globalY = z * 10.0;  // Note: CSV z becomes Y for front/back
        info.globalZ = (surfaceIndex == 4) ? -y * 10.0 : y * 10.0;
        info.rotationAngle = (surfaceIndex == 4) ? 0.0 : 180.0;
    }

    m_leadGlassBlocks[copyNumber] = info;
}

void GeometryManager::RegisterTPCModule(int32_t moduleIndex,
                                         int32_t moduleType,
                                         double x, double y, double z,
                                         double sizeX, double sizeY, double sizeZ,
                                         int32_t driftAxis, int32_t driftSign) {
    TPCModuleInfo info;
    info.moduleIndex = moduleIndex;
    info.moduleType = moduleType;
    info.centerX = x;
    info.centerY = y;
    info.centerZ = z;
    info.sizeX = sizeX;
    info.sizeY = sizeY;
    info.sizeZ = sizeZ;
    info.driftAxis = driftAxis;
    info.driftSign = driftSign;

    // Determine location string
    info.location = (z < 0) ? "front" : "back";

    // Determine position string
    if (moduleType == 2) {  // Type II (horizontal)
        info.position = (y > 0) ? "top" : "bottom";
    } else {  // Type I (vertical)
        std::string vPos = (y > 0) ? "top" : "bottom";
        std::string hPos = (x > 0) ? "right" : "left";
        info.position = vPos + "-" + hPos;
    }

    // Ensure vector is large enough
    if (static_cast<size_t>(moduleIndex) >= m_tpcModules.size()) {
        m_tpcModules.resize(moduleIndex + 1);
    }
    m_tpcModules[moduleIndex] = info;
}

void GeometryManager::RegisterBeampipeSection(int32_t sectionIndex,
                                               const std::string& name,
                                               const std::string& type,
                                               double posZ,
                                               double innerR1, double outerR1,
                                               double innerR2, double outerR2,
                                               double length,
                                               const std::string& material) {
    BeampipeInfo info;
    info.sectionIndex = sectionIndex;
    info.name = name;
    info.type = type;
    info.posZ = posZ;
    info.innerRadius1 = innerR1;
    info.outerRadius1 = outerR1;
    info.innerRadius2 = innerR2;
    info.outerRadius2 = outerR2;
    info.length = length;
    info.material = material;
    info.energyDeposit = 0.0;
    info.nNeutronCaptures = 0;

    m_beampipeSections.push_back(info);
}

// ============================================================================
// Lookup methods
// ============================================================================

const ScintillatorInfo* GeometryManager::GetScintillatorBar(int32_t moduleCopyNumber,
                                                             int32_t barCopyNumber) const {
    auto it = m_scintBarIndex.find({moduleCopyNumber, barCopyNumber});
    if (it != m_scintBarIndex.end()) {
        return &m_scintBars[it->second];
    }
    return nullptr;
}

std::vector<const ScintillatorInfo*> GeometryManager::GetScintillatorBarsOnSurface(int32_t surfaceIndex) const {
    std::vector<const ScintillatorInfo*> result;
    for (const auto& bar : m_scintBars) {
        if (bar.surfaceIndex == surfaceIndex) {
            result.push_back(&bar);
        }
    }
    return result;
}

const LeadGlassInfo* GeometryManager::GetLeadGlassBlock(int32_t copyNumber) const {
    auto it = m_leadGlassBlocks.find(copyNumber);
    if (it != m_leadGlassBlocks.end()) {
        return &it->second;
    }
    return nullptr;
}

std::vector<const LeadGlassInfo*> GeometryManager::GetLeadGlassBlocksOnSurface(int32_t surfaceIndex) const {
    std::vector<const LeadGlassInfo*> result;
    for (const auto& [copyNum, block] : m_leadGlassBlocks) {
        if (block.surfaceIndex == surfaceIndex) {
            result.push_back(&block);
        }
    }
    return result;
}

const TPCModuleInfo* GeometryManager::GetTPCModule(int32_t moduleIndex) const {
    if (moduleIndex >= 0 && static_cast<size_t>(moduleIndex) < m_tpcModules.size()) {
        return &m_tpcModules[moduleIndex];
    }
    return nullptr;
}

const BeampipeInfo* GeometryManager::GetBeampipeSection(int32_t sectionIndex) const {
    for (const auto& sec : m_beampipeSections) {
        if (sec.sectionIndex == sectionIndex) {
            return &sec;
        }
    }
    return nullptr;
}

int32_t GeometryManager::FindTPCModule(double x, double y, double z) const {
    for (const auto& mod : m_tpcModules) {
        // Check if point is within module bounds
        double halfX = mod.sizeX / 2.0;
        double halfY = mod.sizeY / 2.0;
        double halfZ = mod.sizeZ / 2.0;

        if (std::abs(x - mod.centerX) <= halfX &&
            std::abs(y - mod.centerY) <= halfY &&
            std::abs(z - mod.centerZ) <= halfZ) {
            return mod.moduleIndex;
        }
    }
    return -1;  // Not found
}

int32_t GeometryManager::FindLeadGlassBlock(double x, double y, double z, double tolerance) const {
    double minDist = tolerance;
    int32_t bestMatch = -1;

    for (const auto& [copyNum, block] : m_leadGlassBlocks) {
        double dx = x - block.globalX;
        double dy = y - block.globalY;
        double dz = z - block.globalZ;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (dist < minDist) {
            minDist = dist;
            bestMatch = copyNum;
        }
    }
    return bestMatch;
}

std::pair<int32_t, int32_t> GeometryManager::FindScintillatorBar(double x, double y, double z,
                                                                   double tolerance) const {
    double minDist = tolerance;
    std::pair<int32_t, int32_t> bestMatch = {-1, -1};

    for (const auto& bar : m_scintBars) {
        double dx = x - bar.globalX;
        double dy = y - bar.globalY;
        double dz = z - bar.globalZ;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (dist < minDist) {
            minDist = dist;
            bestMatch = {bar.copyNumber, bar.barCopyNumber};
        }
    }
    return bestMatch;
}

// ============================================================================
// Energy deposit recording
// ============================================================================

void GeometryManager::AddLeadGlassEnergyDeposit(int32_t copyNumber, double energy) {
    auto it = m_leadGlassBlocks.find(copyNumber);
    if (it != m_leadGlassBlocks.end()) {
        it->second.energyDeposit += energy;
    }
}

void GeometryManager::AddScintillatorEnergyDeposit(int32_t moduleCopyNumber,
                                                    int32_t barCopyNumber,
                                                    double energy) {
    auto it = m_scintBarIndex.find({moduleCopyNumber, barCopyNumber});
    if (it != m_scintBarIndex.end()) {
        m_scintBars[it->second].energyDeposit += energy;
    }
}

void GeometryManager::AddLeadGlassPhotons(int32_t copyNumber, int32_t nPhotons) {
    auto it = m_leadGlassBlocks.find(copyNumber);
    if (it != m_leadGlassBlocks.end()) {
        it->second.nPhotons += nPhotons;
    }
}

void GeometryManager::AddScintillatorPhotons(int32_t moduleCopyNumber,
                                              int32_t barCopyNumber,
                                              int32_t nPhotons) {
    auto it = m_scintBarIndex.find({moduleCopyNumber, barCopyNumber});
    if (it != m_scintBarIndex.end()) {
        m_scintBars[it->second].nPhotons += nPhotons;
    }
}

void GeometryManager::ClearEventData() {
    for (auto& [copyNum, block] : m_leadGlassBlocks) {
        block.energyDeposit = 0.0;
        block.nPhotons = 0;
        block.nPMTHits = 0;
    }
    for (auto& bar : m_scintBars) {
        bar.energyDeposit = 0.0;
        bar.nPhotons = 0;
    }
}

// ============================================================================
// Statistics and reporting
// ============================================================================

void GeometryManager::PrintStatistics() const {
    std::cout << "\n========================================" << std::endl;
    std::cout << "GeometryManager Statistics" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "TPC Modules:          " << m_tpcModules.size() << std::endl;
    std::cout << "Scintillator Modules: " << m_scintModules.size() << std::endl;
    std::cout << "Scintillator Bars:    " << m_scintBars.size() << std::endl;
    std::cout << "Lead Glass Blocks:    " << m_leadGlassBlocks.size() << std::endl;
    std::cout << "Beampipe Sections:    " << m_beampipeSections.size() << std::endl;
    std::cout << "========================================\n" << std::endl;
}

void GeometryManager::PrintLeadGlassMap() const {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Lead Glass Block Map" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::setw(8) << "CopyNo" << " | "
              << std::setw(8) << "Surface" << " | "
              << std::setw(10) << "X (mm)" << " | "
              << std::setw(10) << "Y (mm)" << " | "
              << std::setw(10) << "Z (mm)" << " | "
              << std::setw(10) << "Edep(MeV)" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // Group by surface
    for (int32_t surf = 0; surf <= 5; ++surf) {
        std::cout << "--- " << Surface::GetName(surf) << " Surface ---" << std::endl;
        auto blocks = GetLeadGlassBlocksOnSurface(surf);
        int count = 0;
        for (const auto* block : blocks) {
            if (count < 5 || block->energyDeposit > 0) {  // Show first 5 and any with energy
                std::cout << std::setw(8) << block->copyNumber << " | "
                          << std::setw(8) << Surface::GetName(block->surfaceIndex) << " | "
                          << std::setw(10) << std::fixed << std::setprecision(1) << block->globalX << " | "
                          << std::setw(10) << block->globalY << " | "
                          << std::setw(10) << block->globalZ << " | "
                          << std::setw(10) << std::setprecision(3) << block->energyDeposit << std::endl;
            }
            count++;
        }
        if (count > 5) {
            std::cout << "  ... (" << (count - 5) << " more blocks on this surface)" << std::endl;
        }
    }
    std::cout << "========================================\n" << std::endl;
}

void GeometryManager::PrintScintillatorMap() const {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Scintillator Module Map" << std::endl;
    std::cout << "========================================" << std::endl;

    for (int32_t surf = 0; surf <= 5; ++surf) {
        auto bars = GetScintillatorBarsOnSurface(surf);
        std::cout << Surface::GetName(surf) << " Surface: "
                  << bars.size() << " bars" << std::endl;
    }
    std::cout << "========================================\n" << std::endl;
}

} // namespace nnbar
