#include "ckalman/c_memory.h"
#include "ckalman/c_math.h"

#include "ccore/c_debug.h"
#include "ccore/c_memory.h"

namespace ncore
{
    namespace nkalman
    {

        nmath::vector_t *memory_t::AllocVector(s32 n)
        {
            u32 const size = sizeof(nmath::vector_t) + sizeof(float) * n;
            ASSERT((u32)(m_memory_current + size - (u8 *)m_memory_base) <= m_memory_size);
            nmath::vector_t *v = (nmath::vector_t *)m_memory_current;
            m_memory_current += sizeof(nmath::vector_t);
            v->m_N    = n;
            v->m_Inc  = 1;
            v->m_data = (float *)m_memory_current;
            m_memory_current += sizeof(float) * n;
            return v;
        }

        nmath::vector_t *memory_t::AllocZeroVector(s32 n)
        {
            nmath::vector_t *v = AllocVector(n);
            nmem::memset(v->m_data, 0, sizeof(float) * n);
            return v;
        }

        nmath::matrix_t *memory_t::AllocMatrix(s32 rows, s32 cols)
        {
            u32 const size = sizeof(nmath::matrix_t) + sizeof(float) * rows * cols;
            ASSERT((u32)(m_memory_current + size - (u8 *)m_memory_base) <= m_memory_size);
            nmath::matrix_t *m = (nmath::matrix_t *)m_memory_current;
            m_memory_current += size;
            m->m_Rows    = rows;
            m->m_Cols    = cols;
            m->m_stride  = cols;
            m->m_capRows = rows;
            m->m_capCols = cols;
            m->m_data    = (float *)m_memory_current;
            return m;
        }

        nmath::matrix_t *memory_t::AllocZeroMatrix(s32 rows, s32 cols)
        {
            nmath::matrix_t *m = AllocMatrix(rows, cols);
            nmem::memset(m->m_data, 0, sizeof(float) * rows * cols);
            return m;
        }

        void memory_t::PushScope()
        {
            ASSERT(m_scope_index + 1 < MEMORY_SCOPE_MAX);
            m_scopes[m_scope_index++] = m_memory_current;
        }

        void memory_t::PopScope()
        {
            ASSERT(m_scope_index - 1 >= 0);
            m_memory_current = m_scopes[--m_scope_index];
#ifdef TARGET_DEBUG
            // In debug mode, clear memory to catch use-after-free bugs
            u32 const size = m_scopes[m_scope_index] - m_memory_current;
            nmem::memset(m_memory_current, 0xCD, (u64)size);
#endif
        }

        void memory_t::Reset() { m_memory_current = (u8 *)m_memory_base; }

    }  // namespace nkalman
}  // namespace ncore
