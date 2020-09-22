#if !defined (INCLUDE_PHRQALLOC_H)
#define INCLUDE_PHRQALLOC_H

#if defined (USE_PHRQ_ALLOC)

#if !defined(NDEBUG)
void *PHRQ_malloc(size_t, const char *, int);
void *PHRQ_calloc(size_t, size_t, const char *, int);
void *PHRQ_realloc(void *, size_t, const char *, int);
#else
extern void *PHRQ_malloc(size_t);
extern void *PHRQ_calloc(size_t, size_t);
extern void *PHRQ_realloc(void *, size_t);
#endif

extern void PHRQ_free(void *);
extern void PHRQ_free_all(void);

#if !defined(NDEBUG)
#define   PHRQ_malloc(s)         PHRQ_malloc(s, __FILE__, __LINE__)
#define   PHRQ_calloc(c, s)      PHRQ_calloc(c, s, __FILE__, __LINE__)
#define   PHRQ_realloc(p, s)     PHRQ_realloc(p, s, __FILE__, __LINE__)
#endif

#else /* defined (USE_PHRQ_ALLOC) */

#if !defined(NDEBUG)
void *PHRQ_malloc(size_t, const char *, int);
void *PHRQ_calloc(size_t, size_t, const char *, int);
void *PHRQ_realloc(void *, size_t, const char *, int);
#else
extern void *PHRQ_malloc(size_t);
extern void *PHRQ_calloc(size_t, size_t);
extern void *PHRQ_realloc(void *, size_t);
#endif
void PHRQ_free(void *);
void PHRQ_free_all(void);

#if !defined(NDEBUG)
#define   PHRQ_malloc(s)         PHRQ_malloc(s, __FILE__, __LINE__)
#define   PHRQ_calloc(c, s)      PHRQ_calloc(c, s, __FILE__, __LINE__)
#define   PHRQ_realloc(p, s)     PHRQ_realloc(p, s, __FILE__, __LINE__)
#endif

#endif /* defined (USE_PHRQ_ALLOC) */

#endif /* !defined (INCLUDE_PHRQALLOC_H) */
