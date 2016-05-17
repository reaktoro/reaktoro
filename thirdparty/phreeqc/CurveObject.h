#if !defined(CURVEOBJECT_H_INCLUDED)
#define CURVEOBJECT_H_INCLUDED
#include <vector>
#include <string>
#include "phrqtype.h"
class CurveObject
{

public:
	CurveObject();
	~CurveObject();
	void Set_id(std::string s)
	{
		this->id = s;
	}
	std::string &Get_id(void)
	{
		return this->id;
	}
	void Set_color(std::string s)
	{
		this->color = s;
	}
	std::string &Get_color(void)
	{
		return this->color;
	}
	void Set_symbol(std::string s)
	{
		this->symbol = s;
	}
	std::string &Get_symbol(void)
	{
		return this->symbol;
	}
	void Set_symbol_size(LDBLE f)
	{
		this->symbol_size = f;
	}
	LDBLE Get_symbol_size(void)
	{
		return this->symbol_size;
	}
	void Set_line_w(LDBLE f)
	{
		this->line_w = f;
	}
	LDBLE Get_line_w(void)
	{
		return this->line_w;
	}
	void Set_y_axis(int f)
	{
		this->y_axis = f;
	}
	std::vector<LDBLE> & Get_x()
	{
		return this->x;
	}
	std::vector<LDBLE> & Get_y()
	{
		return this->y;
	}
	int Get_y_axis()
	{
		return this->y_axis;
	}

protected:
	std::vector<LDBLE> x, y;
	std::string id, color, symbol;
	int y_axis; 
	LDBLE line_w, symbol_size;

public:

};

#endif // !defined(CURVEOBJECT_H_INCLUDED)
