 #pragma once
#include "ThirdPartyHeadersBegin.h"
#  include <string>
#include "ThirdPartyHeadersEnd.h"
#include "IJK.h"
#include "fileStuff.h"
#include "SzlFileLoader.h"
namespace tecplot { namespace ___3933 { class ___1388 { public: static inline bool validIsAscii(___372 ___2002) { return VALID_BOOLEAN(___2002); }; static inline bool validDataFileType(DataFileType_e ___844) { return VALID_ENUM(___844, DataFileType_e); } static inline bool validIJKSubzoneSize(___1844 const& ijkSubzoneSize) { return ijkSubzoneSize.___2067() && ijkSubzoneSize.blockSize()<=___2090::MAX_ITEM_OFFSET+1; }; static inline bool validFESubzoneSize(___2090::ItemOffset_t feSubzoneSize) { return feSubzoneSize>0 && feSubzoneSize<=___2090::MAX_ITEM_OFFSET+1; }; static inline bool validFileVersion(uint32_t fileVersion) { return fileVersion >= SZPLT_MIN_READ_VERSION; } private: ___372                 m_isAscii; DataFileType_e            m_dataFileType; ___1844                       m_maxIJKSubzoneSize; ___2090::ItemOffset_t m_maxFESubzoneSize; uint32_t                  m_fileVersion; uint32_t                  m_codeRevision; public: ___1388() { invalidate(); } ___1388( ___372 const                 ___2002, DataFileType_e                  ___844, ___1844 const&                      ijkSubzoneSize, ___2090::ItemOffset_t const feSubzoneSize) { invalidate(); ___3494(___2002); setDataFileType(___844); setMaxIJKSubzoneSize(ijkSubzoneSize); setMaxFESubzoneSize(feSubzoneSize); setFileVersion(SZPLT_CUR_WRITE_VERSION); } inline void invalidate(void) { m_isAscii = ___372(-1); m_dataFileType = ___847; m_maxIJKSubzoneSize.invalidate(); m_maxFESubzoneSize = ___2090::INVALID_ITEM_OFFSET; m_fileVersion = 0; m_codeRevision = 0; }
 #if 0
inline bool operator ==(___1388 const& ___1308) const { return ___2002() == ___1308.___2002() && ___844() == ___1308.m_dataFileType() && ___1757() == ___1308.___1757() && ___1756() == ___1308.___1756() && getFileVersion() == ___1308.getFileVersion() && getRevision() == ___1308.getRevision(); } inline bool operator !=(___1388 const& ___1308) const { return !(*this == ___1308); }
 #endif
inline void ___3494(___372 ___2002) { REQUIRE(validIsAscii(___2002)); m_isAscii = ___2002; } inline void setDataFileType(DataFileType_e ___844) { REQUIRE(validDataFileType(___844)); m_dataFileType = ___844; } inline void setMaxIJKSubzoneSize(___1844 const& maxIJKSubzoneSize) { REQUIRE(validIJKSubzoneSize(maxIJKSubzoneSize)); m_maxIJKSubzoneSize = maxIJKSubzoneSize; } inline void setMaxFESubzoneSize(___2090::ItemOffset_t maxFESubzoneSize) { REQUIRE(validFESubzoneSize(maxFESubzoneSize)); m_maxFESubzoneSize = maxFESubzoneSize; } inline void setFileVersion(uint32_t fileVersion) { m_fileVersion = fileVersion; } inline void setCodeRevision(uint32_t codeRevision) { m_codeRevision = codeRevision; } inline ___372 ___2002() const { ENSURE(VALID_BOOLEAN(m_isAscii)); return m_isAscii; } inline DataFileType_e ___844() const { ENSURE(VALID_ENUM(m_dataFileType, DataFileType_e)); return m_dataFileType; } inline ___1844 const& ___1757() const { ENSURE(validIJKSubzoneSize(m_maxIJKSubzoneSize)); return m_maxIJKSubzoneSize; } inline ___2090::ItemOffset_t ___1756() const { ENSURE(validFESubzoneSize(m_maxFESubzoneSize)); return m_maxFESubzoneSize; } inline uint32_t getFileVersion() const { ENSURE(validFileVersion(m_fileVersion)); return m_fileVersion; } inline uint32_t getCodeRevision() const { return m_codeRevision; } }; }}
