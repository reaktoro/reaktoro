#define INCLUDE_PHRQALLOC_H
#include "Phreeqc.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if defined(PHREEQCI_GUI)
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif
#if defined(USE_PHRQ_ALLOC)
/* ---------------------------------------------------------------------- */
#if !defined(NDEBUG)
void * Phreeqc::
PHRQ_malloc(size_t size, const char *szFileName, int nLine)
#else
void *Phreeqc::
PHRQ_malloc(size_t size)
#endif
/* ---------------------------------------------------------------------- */
{
	PHRQMemHeader *p;

	assert((s_pTail == NULL) || (s_pTail->pNext == NULL));

	p = (PHRQMemHeader *) malloc(sizeof(PHRQMemHeader) + size);

	if (p == NULL)
		return NULL;
#if !defined(NDEBUG)
	memset(p, 0, sizeof(PHRQMemHeader) + size);
#endif
	p->pNext = NULL;

	if ((p->pPrev = s_pTail) != NULL)
	{
		s_pTail->pNext = p;
	}

	p->size = sizeof(PHRQMemHeader) + size;
#if !defined(NDEBUG)
	p->szFileName = (char *) malloc(strlen(szFileName) + 1);
	if (p->szFileName)
		strcpy(p->szFileName, szFileName);
	p->nLine = nLine;
#endif

	s_pTail = p;
	p++;
	return ((void *) (p));
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
PHRQ_free(void *ptr)
/* ---------------------------------------------------------------------- */
{
	PHRQMemHeader *p;

	assert((s_pTail == NULL) || (s_pTail->pNext == NULL));

	if (ptr == NULL)
		return;

	p = (PHRQMemHeader *) ptr - 1;

	if (p->pNext != NULL)
	{
		p->pNext->pPrev = p->pPrev;
	}
	else
	{
		/* Handle special case when (p == s_pTail) */
		assert(s_pTail != NULL);
		assert(p == s_pTail);
		s_pTail = p->pPrev;
	}

	if (p->pPrev)
	{
		p->pPrev->pNext = p->pNext;
	}

#if !defined(NDEBUG)
	free(p->szFileName);
#endif

	free(p);
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
PHRQ_free_all(void)
/* ---------------------------------------------------------------------- */
{
	assert((s_pTail == NULL) || (s_pTail->pNext == NULL));

	std::ostringstream ostrm;

	if (s_pTail == NULL)
	{
#if !defined(NDEBUG)
		output_msg("No memory leaks\n");
#endif
		return;
	}
	while (s_pTail->pPrev != NULL)
	{
		s_pTail = s_pTail->pPrev;
#if !defined(NDEBUG)
		ostrm.clear();
		ostrm << s_pTail->pNext->szFileName << "(" << s_pTail->pNext->nLine;
		ostrm << ") " << (void *) (s_pTail->pNext + 1) << ": freed in PHRQ_free_all\n";
		output_msg(ostrm.str().c_str());
		free(s_pTail->pNext->szFileName);
#endif
		free(s_pTail->pNext);
	}

#if !defined(NDEBUG)
	ostrm.clear();
	ostrm <<  s_pTail->szFileName << "(" << s_pTail->nLine;
	ostrm << ") " << (void *) (s_pTail + 1) << ": freed in PHRQ_free_all\n";
	output_msg(ostrm.str().c_str());
	free(s_pTail->szFileName);
#endif
	free(s_pTail);
	s_pTail = NULL;
}

/* ---------------------------------------------------------------------- */
void * Phreeqc::
PHRQ_calloc(size_t num, size_t size
#if !defined(NDEBUG)
			, const char *szFileName, int nLine
#endif
	)
/* ---------------------------------------------------------------------- */
{
	PHRQMemHeader *p;

	assert((s_pTail == NULL) || (s_pTail->pNext == NULL));

	p = (PHRQMemHeader *) malloc(sizeof(PHRQMemHeader) + size * num);

	if (p == NULL)
		return NULL;

	p->pNext = NULL;

	if ((p->pPrev = s_pTail) != NULL)
	{
		s_pTail->pNext = p;
	}

	p->size = sizeof(PHRQMemHeader) + size * num;

#if !defined(NDEBUG)
	p->szFileName = (char *) malloc(strlen(szFileName) + 1);
	if (p->szFileName)
		strcpy(p->szFileName, szFileName);
	p->nLine = nLine;
#endif

	s_pTail = p;
	p++;
	return memset(p, 0, size * num);
}

/* ---------------------------------------------------------------------- */
void * Phreeqc::
PHRQ_realloc(void *ptr, size_t size
#if !defined(NDEBUG)
			 , const char *szFileName, int nLine
#endif
	)
/* ---------------------------------------------------------------------- */
{
	PHRQMemHeader *p;
	size_t new_size;
	size_t old_size;

	if (ptr == NULL)
	{
		return PHRQ_malloc(size
#if !defined(NDEBUG)
						   , szFileName, nLine
#endif
			);
	}

	assert((s_pTail == NULL) || (s_pTail->pNext == NULL));

	p = (PHRQMemHeader *) ptr - 1;

	new_size = sizeof(PHRQMemHeader) + size;

	old_size = p->size;
	p = (PHRQMemHeader *) realloc(p, new_size);
	if (p != NULL)
	{
		p->size = new_size;
#if !defined(NDEBUG)
		if (new_size > old_size)
		{
			memset((char *) p + old_size, 0, new_size - old_size);
		}
#endif
	}

	if (p == NULL)
		return NULL;

	if (p->pPrev != NULL)
	{
		p->pPrev->pNext = p;
	}

	if (p->pNext != NULL)
	{
		p->pNext->pPrev = p;
	}
	else
	{
		s_pTail = p;
	}

#if !defined(NDEBUG)
	free(p->szFileName);
	p->szFileName = (char *) malloc(strlen(szFileName) + 1);
	if (p->szFileName)
		strcpy(p->szFileName, szFileName);
	p->nLine = nLine;
#endif

	p++;
	return ((void *) (p));
}
#else /* USE_PHRQ_ALLOC */
/* ---------------------------------------------------------------------- */
void *Phreeqc::
#if !defined(NDEBUG)
PHRQ_malloc(size_t size, const char *szFileName, int nLine)
#else
PHRQ_malloc(size_t size)
#endif
/* ---------------------------------------------------------------------- */
{
#if !defined(NDEBUG) && defined(WIN32_MEMORY_DEBUG)
	return _malloc_dbg(size, _NORMAL_BLOCK, szFileName, nLine);
#else
	return malloc(size);
#endif
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
PHRQ_free(void *ptr)
/* ---------------------------------------------------------------------- */
{
#if !defined(NDEBUG) && defined(WIN32_MEMORY_DEBUG)
	_free_dbg(ptr, _NORMAL_BLOCK);
#else
	free(ptr);
#endif
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
PHRQ_free_all(void)
/* ---------------------------------------------------------------------- */
{
}

/* ---------------------------------------------------------------------- */
void * Phreeqc::
#if !defined(NDEBUG)
PHRQ_calloc(size_t num, size_t size, const char *szFileName, int nLine)
#else
PHRQ_calloc(size_t num, size_t size)
#endif
/* ---------------------------------------------------------------------- */
{
#if !defined(NDEBUG) && defined(WIN32_MEMORY_DEBUG)
	return _calloc_dbg(num, size, _NORMAL_BLOCK, szFileName, nLine);
#else
	return calloc(num, size);
#endif
}

/* ---------------------------------------------------------------------- */
void * Phreeqc::
#if !defined(NDEBUG)
PHRQ_realloc(void *ptr, size_t size, const char *szFileName, int nLine)
#else
PHRQ_realloc(void *ptr, size_t size)
#endif
/* ---------------------------------------------------------------------- */
{
#if !defined(NDEBUG) && defined(WIN32_MEMORY_DEBUG)
	return _realloc_dbg(ptr, size, _NORMAL_BLOCK, szFileName, nLine);
#else
	return realloc(ptr, size);
#endif
}
#endif /* USE_PHRQ_ALLOC */
