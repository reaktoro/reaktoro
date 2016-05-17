#if !defined(CHARTOBJECT_H_INCLUDED)
#define CHARTOBJECT_H_INCLUDED
#if defined MULTICHART
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include "CurveObject.h"
#include "NumKeyword.h"

#include <float.h>
class Phreeqc;

class ChartObject:public cxxNumKeyword
{

  public:
	ChartObject(PHRQ_io *io=NULL);
	ChartObject(int i, PHRQ_io *io=NULL);
	~ChartObject();

	enum chart_batch_type
	{ ChO_NO_BATCH   = -1, 
	  ChO_BATCH_ONLY = 0,
	  ChO_EMF        = 1,
	  ChO_PNG        = 2, 
	  ChO_JPG        = 3,
	  ChO_GIF        = 4,
	  ChO_TIFF       = 5, 
	  ChO_BMP        = 6
	};

	bool Get_new_ug()const
	{
		return this->new_ug;
	}
	void Set_new_ug(bool b)
	{
		this->new_ug = b;
	}
	void Set_FirstCallToUSER_GRAPH(bool b)
	{
		this->FirstCallToUSER_GRAPH = b;
	}
	bool Get_FirstCallToUSER_GRAPH()const
	{
		return this->FirstCallToUSER_GRAPH;
	}
	int Get_update_time_chart()const
	{
		return (this->update_time_chart);
	}
	int Get_PanelHeight()const
	{
		return (this->PanelHeight);
	}
	int Get_PanelWidth()const
	{
		return (this->PanelWidth);
	}
	std::string &Get_chart_title()
	{
		return this->chart_title;
	}
	const std::string &Get_chart_title()const
	{
		return this->chart_title;
	}
	std::string Get_batch_fn()
	{
		return this->batch_fn;
	}
	void Set_batch_fn(std::string fn)
	{
		this->batch_fn = fn;
	}
	std::vector<std::string> &Get_axis_titles()
	{
		return this->axis_titles;
	}
	const std::vector<std::string> &Get_axis_titles()const
	{
		return this->axis_titles;
	}
	LDBLE *Get_axis_scale_x()
	{
		return this->axis_scale_x;
	}
	const LDBLE *Get_axis_scale_x()const
	{
		return this->axis_scale_x;
	}
	LDBLE *Get_axis_scale_y()
	{
		return this->axis_scale_y;
	}
	const LDBLE *Get_axis_scale_y()const
	{
		return this->axis_scale_y;
	}
	LDBLE *Get_axis_scale_y2()
	{
		return this->axis_scale_y2;
	}
	const LDBLE *Get_axis_scale_y2()const
	{
		return this->axis_scale_y2;
	}
	int Get_chart_type()const
	{
		return this->chart_type;
	}
	bool Get_graph_initial_solutions()const
	{
		return this->graph_initial_solutions;
	}
	void Set_graph_initial_solutions(bool val)
	{
		this->graph_initial_solutions = val;
	}
	bool Get_connect_simulations()const
	{
		return this->connect_simulations;
	}
	void Set_connect_simulations(bool val)
	{
		this->connect_simulations = val;
	}
	void Set_colnr(int i)
	{
		this->colnr = i;
	}
	int Get_colnr()const
	{
		return (this->colnr);
	}
	void Set_ColumnOffset(int i)
	{
		this->ColumnOffset = i;
	}
	int Get_ColumnOffset()const
	{
		return (this->ColumnOffset);
	}	
	void Set_AddSeries(bool b)
	{
		this->AddSeries = b;
	}
	bool Get_AddSeries()const
	{
		return this->AddSeries;
	}
	std::vector< std::string > Get_csv_file_names()const
	{
		return this->csv_file_names;
	}
	void Get_csv_file_names(std::vector< std::string > names)
	{
		this->csv_file_names = names;
	}
	void Set_prev_advection_step(int i)
	{
		this->prev_advection_step = i;
	}
	int Get_prev_advection_step()const
	{
		return (this->prev_advection_step);
	}
	void Set_prev_transport_step(int i)
	{
		this->prev_transport_step = i;
	}
	int Get_prev_transport_step()const
	{
		return (this->prev_transport_step);
	}
	void Set_prev_sim_no(int i)
	{
		this->prev_sim_no = i;
	}
	int Get_prev_sim_no(void)const
	{
		return this->prev_sim_no;
	}
	void Set_end_timer(bool b)
	{
		this->end_timer = b;
	}
	bool Get_end_timer()const
	{
		return this->end_timer;
	}
	void Set_done(bool b)
	{
		this->done = b;
	}
	bool Get_done()const
	{
		return this->done;
	}
	std::vector<CurveObject *> &Get_CurvesCSV()
	{
		return this->CurvesCSV;
	}
	const std::vector<CurveObject *> &Get_CurvesCSV()const
	{
		return this->CurvesCSV;
	}
	std::vector<CurveObject *> &Get_Curves()
	{
		return this->Curves;
	}
	const std::vector<CurveObject *> &Get_Curves()const
	{
		return this->Curves;
	}
	void Set_curve_added(bool tf)
	{
		this->curve_added = tf;
	}
	bool Get_curve_added()const
	{
		return this->curve_added;
	}
	void Set_point_added(bool tf)
	{
		this->point_added = tf;
	}
	bool Get_point_added()const
	{
		return this->point_added;
	}
	struct rate *Get_user_graph()
	{
		return this->user_graph;
	}
	const struct rate *Get_user_graph()const
	{
		return this->user_graph;
	}
	std::list<std::string> &Get_rate_command_list()
	{
		return this->rate_command_list;
	}
	const std::list<std::string> &Get_rate_command_list()const
	{
		return this->rate_command_list;
	}
	void Set_rate_new_def(bool tf);

	bool Get_rate_new_def()const
	{
		return this->rate_new_def;
	}
	void Set_graph_x(LDBLE d)
	{
		this->graph_x = d;
	}
	LDBLE Get_graph_x()const
	{
		return this->graph_x;
	}
	std::map<int, LDBLE> &Get_graph_y()
	{
		return this->graph_y;
	}
	const std::map<int, LDBLE> &Get_graph_y()const
	{
		return this->graph_y;
	}
	std::map<int, bool> &Get_secondary_y()
	{
		return this->secondary_y;
	}
	const std::map<int, bool> &Get_secondary_y()const
	{
		return this->secondary_y;
	}
	std::vector<CurveObject> &Get_new_plotxy_curves()
	{
		return this->new_plotxy_curves;
	}
	const std::vector<CurveObject> &Get_new_plotxy_curves()const
	{
		return this->new_plotxy_curves;
	}
	std::vector<std::string> &Get_new_headings()
	{
		return this->new_headings;
	}
	const std::vector<std::string> &Get_new_headings()const
	{
		return this->new_headings;
	}
	void Set_active(bool tf)
	{
		this->active = tf;
	}
	bool Get_active()const
	{
		return this->active;
	}
	void Set_detach(bool tf)
	{
		this->detach = tf;
	}
	bool Get_detach()const
	{
		return this->detach;
	}
	bool Get_form_started()const
	{
		return this->form_started;
	}
	void Set_phreeqc(Phreeqc * ptr)
	{
		this->phreeqc_ptr = ptr;
	}
	Phreeqc * Get_phreeqc()
	{
		return this->phreeqc_ptr;
	}
	const Phreeqc * Get_phreeqc()const
	{
		return this->phreeqc_ptr;
	}
	const std::list<std::string>& Get_rate_command_list_original()const
	{
		return this->rate_command_list_original;
	}
	bool Set_axis_scale(std::vector<std::string>, std::vector<int> types, std::ostringstream &);
	bool Set_axis_scale(CParser & parser);
	bool Read(CParser & parser);
	bool OpenCSVFile(std::string fn);
	static CurveObject ExtractCurveInfo(std::string & str_line);
	void Set_rate_struct(void);
	void PlotXY(std::string x, std::string y);
    bool start_chart(void);
	ZedGraph::SymbolType Return_SymbolType(std::string);
	void SaveCurvesToFile(std::string &);
	void Rate_free(void);
	void Initialize_graph_pts(void);
	void Finalize_graph_pts(void);
	void Get_legal_symbol(std::string &sym);
	void Get_legal_symbol_csv(std::string &sym);
	void Get_color_string(std::string &color);
	void Get_color_string_csv(std::string &color);
	void Add_new_series(void);
	void Add_curve(bool plotxy, std::string id = "", 
					   LDBLE line_width = 1.0, 
					   std::string symbol = "",
					   LDBLE symbol_size = 6.0, 
					   int y_axis = 1,
					   std::string color = "");
	void dump(std::ostream & s_oss, unsigned int indent);

	chart_batch_type Get_batch(void) {return this->batch;}
	void Set_batch(chart_batch_type bt) {this->batch = bt;}
	bool Get_batch_background(void) {return this->batch_background;}
	void Set_batch_background(bool tf) {this->batch_background = tf;}

	bool Get_batch_grid(void) {return this->batch_grid;}
	void Set_batch_grid(bool tf) {this->batch_grid = tf;}
  protected:

	bool new_ug;
	bool FirstCallToUSER_GRAPH;

	int update_time_chart;			/* milliseconds, maybe read */
	int PanelHeight;
	int PanelWidth;
	std::map<std::string, int> Symbol_map;
	std::vector<std::string> Color_vector; 
	std::string chart_title;
	std::vector<std::string> axis_titles;
	LDBLE axis_scale_x[5];
	LDBLE axis_scale_y[5];
	LDBLE axis_scale_y2[5];

	int chart_type;
	bool graph_initial_solutions;
	bool connect_simulations;
	int shifts_as_points;
	int colnr;
	int ColumnOffset;
	bool AddSeries;
	
	int prev_advection_step;
	int prev_transport_step;

	int prev_sim_no;

	bool end_timer;
	bool done;

	std::vector<std::string> csv_file_names;
	std::vector<CurveObject *> CurvesCSV;
	std::vector<CurveObject *> Curves;
	bool curve_added;
	bool point_added;
	
	struct rate *user_graph;
	// C++ for rate struct
	std::string rate_name;
	std::list<std::string> rate_command_list;
	std::list<std::string> rate_command_list_original;
	bool rate_new_def;

	int default_symbol;
	int default_symbol_csv;
	int default_color;
	int default_color_csv;

	// temporary storage before stored graph_x/y/sy data are stored in curves
	// Initialize_graph_pts and Finalize_graph_pts use this storage.
	LDBLE graph_x;
	std::map<int, LDBLE> graph_y;
	std::map<int, bool> secondary_y;

	// temporary plotxy curve definitions before stored in curves
	// a plotxy curve is copied to Curves when cmdplotxy is encountered
	// this keeps order correct between plotxy and graph_x/y/sy
	std::vector<CurveObject> new_plotxy_curves;

	// temporary headings until stored during basic_run
	std::vector<std::string> new_headings;
	std::vector<std::string> headings_original;
	bool active;
	bool detach;
	bool form_started;
	class Phreeqc * phreeqc_ptr;

	bool batch_background;
	bool batch_grid;
	std::string batch_fn;
	chart_batch_type batch;
	const static std::vector < std::string > vopts;

  public:
	int usingResource;

};
#endif // MULTICHART
#endif // !defined(CHARTOBJECT_H_INCLUDED)
