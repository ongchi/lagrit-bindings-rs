#include <stdlib.h>

// File control functions
extern void fc_fflush_and_sync(const char *file);
extern void fc_fclose(const char *file);

// Memory management functions
extern void fc_mmfindbk(size_t *byte_len, size_t *cell_len, const char *aname,
                        const char *pname, int8_t *aptr, int32_t *arr_len,
                        int32_t *status);
extern void fc_mmrelprt(const char *pname, int32_t *status);

// Lagrit functions
extern void fc_initlagrit(int32_t *mode, const char *log_file,
                          const char *batch_file);
extern void fc_dotask(const char *cmd, int32_t *status);
extern void fc_attr_len(const char *aname, const char *pname, int32_t *arr_len,
                        int32_t *status);
extern void fc_cmo_get_index(const char *cmo_name, int32_t *idx,
                             int32_t *status);
extern void fc_cmo_get_name(const char *cmo_name, int32_t *status);
extern void fc_cmo_get_mesh_type(const char *cmo_name, const char *mesh_type,
                                 int32_t *imesh_type, int32_t *status);
