#pragma once
#include <Windows.h>
#include <cassert>	

namespace zdg_ui2 {
	using namespace System;
	//using namespace System::ComponentModel;
	using namespace System::Resources;
	using namespace System::Windows::Forms;
	using namespace System::Drawing;
	using namespace System::Drawing::Imaging;
	using namespace System::Threading;
	using namespace ZedGraph;
							//using namespace System::Runtime::InteropServices;
							//[DllImport("gdi32.dll")]
							////extern IntPtr CopyEnhMetaFileA(IntPtr hemfSrc, System::Text::StringBuilder* hNULL);
							//extern IntPtr CopyEnhMetaFileA(IntPtr hemfSrc, String^ hNULL);
// Form1 is only used with MULTICHART
	public ref class ChartObj : public System::Object
	{
	public: ChartObject* chartobject_ptr;
	public:	ChartObj(ChartObject* ptr)
			{
				this->chartobject_ptr = ptr;
			}
	};

	public ref class Form1  : public System::Windows::Forms::Form
	{
	public:	long int tickStart;
	public:	Form1(ChartObject *ptr)
			{
				this->chartobject_ptr = ptr;
				if (Phreeqc* ptr = this->chartobject_ptr->Get_phreeqc())
				{
					ptr->Get_chart_handler().Increment_active_charts();
				}
				InitializeComponent();
				col_use = 0;
				symbol_use = 0;
				Y2 = false;
				phreeqc_done = false;
				Y2show = false;	
				background = true;
				hints = true;
				grid = true;
			}
			static void ThreadForm(Object^ data)
			{
				ChartObject *ptr = ((ChartObj^)(data))->chartobject_ptr;
				Form1 ^myForm = gcnew Form1(ptr);
				myForm->ShowDialog();
				myForm->~Form1();
			}
	private: bool phreeqc_done;
			 
	private: void SetSize()
			 {
				 zg1->Location = Point( 0, 0 );
				 // Leave a small margin around the outside of the control
				 zg1->Size = System::Drawing::Size( ClientRectangle.Width - 0,
					 ClientRectangle.Height - 0 );
			 }

			System::Void MyFormClosingEventHandler(
					System::Object^ sender, 
				System::Windows::Forms::FormClosingEventArgs ^e)
			{
				ChartObject *chart = this->chartobject_ptr;
				if (chart != NULL) 
				{
					chart->Set_done(true);
					System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
					if (Phreeqc* ptr = chart->Get_phreeqc())
					{
						ptr->Get_chart_handler().Decrement_active_charts();
						//chart->Set_phreeqc(NULL);
					}
				}
				this->chartobject_ptr = NULL;
			}

			 System::Void Form1_Load(System::Object ^sender, System::EventArgs ^e)
			 {
				 CreateGraph( zg1 );
				 SetSize();
			 }

			 System::Void Form1_Resize(System::Object ^sender, System::EventArgs ^e)
			 {
				 SetSize();
			 }

			 //static bool LogX, LogY, LogY2;
			 bool LogX, LogY, LogY2;
	private: bool check_neg_log( int i, int i2)
			 {
				 ChartObject *chart = this->chartobject_ptr;
				 if (chart == NULL) return false;
				 std::vector<CurveObject *> &Curves = chart->Get_Curves();
				 if (LogX && chart->Get_axis_scale_x()[4] == 10.0 && 
					 Curves[i]->Get_x()[i2] <= 0)
				 {
					 if (Phreeqc* ptr = chart->Get_phreeqc())
					 {
						 ptr->warning_msg("Obtained x_value <= 0, removing point...");
					 }
					 //axis_scale_x[4] = NA; /* if reverting to linear... */
					 //LogX = false;
					 return true;
				 }
				 if (Curves[i]->Get_y()[i2] <= 0 && 
					 (chart->Get_axis_scale_y()[4] == 10.0 || 
					 chart->Get_axis_scale_y2()[4] == 10.0))
				 {
					 if (Curves[i]->Get_y_axis() == 2 && LogY2)
					 {
						 if (Phreeqc* ptr = chart->Get_phreeqc())
						 {
							 ptr->warning_msg("Obtained sy_value <= 0, removing point......");
						 }
						 //axis_scale_y2[4] = NA;
						 //LogY2 = false;
						 return true;
					 }
					 else if (LogY)
					 {
						 if (Phreeqc* ptr = chart->Get_phreeqc())
						 {
							 ptr->warning_msg("Obtained y_value <= 0, removing point......");
						 }
						 //axis_scale_y[4] = NA;
						 //LogY = false;
						 return true;
					 }
				 }
				 return false;
			 }

	private: PointPairList ^list;
			 int col_use, symbol_use;
			 bool Y2, Y2show;
	//		 static cli::array<String^> ^ColorList = {"Red", "Green", "Blue", "Orange", "Magenta", "Yellow", "Black", "Cyan", "Brown", "Lime", "Gray" };
			 static cli::array<String^> ^ColorList = {"Red", "Green", "Blue", "Orange", "Magenta", "Black", "Cyan", "Brown", "Lime", "Gray" };
			 bool background, hints, grid;

			 ZedGraph::GraphObjList ^GOL_no_hints;
			 ZedGraph::GraphObjList ^GOL_hints;

			 void DefineCurves(GraphPane ^myPane, int init)
			 {
				 ChartObject *chart = this->chartobject_ptr;
				 if (chart == NULL) 
				 {
					 return;
				 }
 				 std::vector<CurveObject *> Curves; 
				 size_t i;
				 for (i = 0; i < chart->Get_CurvesCSV().size(); i++)
				 {
					 Curves.push_back(chart->Get_CurvesCSV()[i]);
				 }
				 for (i = 0; i < chart->Get_Curves().size(); i++)
				 {
					 Curves.push_back(chart->Get_Curves()[i]);
				 }

				 chart->Set_curve_added(false);

				 // Set the titles and axis labels
				 myPane->Title->Text = gcnew String(chart->Get_chart_title().c_str());
				 if (chart->Get_axis_titles().size() > 0)
					 myPane->XAxis->Title->Text = gcnew String(chart->Get_axis_titles()[0].c_str());
				 if (chart->Get_axis_titles().size() > 1)
					 myPane->YAxis->Title->Text = gcnew String(chart->Get_axis_titles()[1].c_str());
				 if (chart->Get_axis_titles().size() > 2)
					 myPane->Y2Axis->Title->Text = gcnew String(chart->Get_axis_titles()[2].c_str());

				 LineItem ^myCurve;

				 Color col;

				 String ^s_t;
				 if (chart->Get_axis_scale_x()[4] == 10.0) LogX = true;
				 else LogX = false;
				 if (chart->Get_axis_scale_y()[4] == 10.0) LogY = true;
				 else LogY = false;
				 if (chart->Get_axis_scale_y2()[4] == 10.0) LogY2 = true;
				 else LogY2 = false;

				 //Rewrite all curves
				 zg1->GraphPane->CurveList->Clear();
				 for (size_t i = 0; i < Curves.size(); i++)
				 {
					 // even curves with no data
					 //if (Curves[i]->Get_x().size() == 0) continue;
					 list = gcnew PointPairList();
					 if (Curves[i]->Get_y_axis() == 2)
					 {
						 Y2 = true;
						 Y2show = true;
					 }
					 else
						 Y2 = false;
					 for (int i2 = 0; (i2 < (int) Curves[i]->Get_x().size()); i2++)
					 {
						 if ((LogX && Curves[i]->Get_x()[i2] <=0)
							 || (LogY && !Y2 && Curves[i]->Get_y()[i2] <=0)
							 || (LogY2 && Y2 && Curves[i]->Get_y()[i2] <=0))
							 continue;
						 else
							 list->Add( Curves[i]->Get_x()[i2], Curves[i]->Get_y()[i2] );
					 }

					 col = Color::FromName(gcnew String(Curves[i]->Get_color().c_str()));
					 if (!col.IsKnownColor)
					 {
						 col = Color::FromName(ColorList[col_use]);
						 std::string newcol;
						 ToString(col.ToString(), newcol);
						 Utilities::replace("Color [","",newcol);
						 Utilities::replace("]","",newcol);
						 Curves[i]->Set_color(newcol);

					 }
					 if (++col_use > 6) col_use = 0;

					 SymbolType symb = chart->Return_SymbolType
						 (Curves[i]->Get_symbol());

					 // id
					 s_t = gcnew String(Curves[i]->Get_id().c_str());


					 // Add curve to chart
					 myCurve = myPane->AddCurve( s_t, list, col, symb );


					 // Curve with no points is invisible
					 if (Curves[i]->Get_x().size() == 0) 
					 {
						 myCurve->IsVisible = false;
						 myCurve->Label->IsVisible = false;
					 }

					 if (Curves[i]->Get_line_w() > 0.0)
						 myCurve->Line->Width = (float) Curves[i]->Get_line_w();
					 else
						 myCurve->Line->IsVisible = false;
					 /* hmm... dash/dot don`t display well */
					 // myCurve->Line->Style = System::Drawing::Drawing2D::DashStyle::Dot;
					 myCurve->Symbol->Fill = gcnew Fill( Color::FromName("White") );
					 if (Curves[i]->Get_symbol_size() > 0.0)
						 myCurve->Symbol->Size = (float) Curves[i]->Get_symbol_size();
					 else
						 myCurve->Symbol->IsVisible = false;
					 myCurve->Symbol->Border->Width = (float) Curves[i]->Get_line_w();
					 if (Y2)
						 myCurve->IsY2Axis = true;
					 delete list;
				 }

				 if (Y2show)
					 myPane->Legend->Position = ZedGraph::LegendPos::TopCenter;
				 else
					 myPane->Legend->Position = ZedGraph::LegendPos::Right;
				 myPane->Legend->FontSpec->Size = 12;
				 myPane->Legend->FontSpec->IsBold = false;

				 // Show the x axis grid
				 myPane->XAxis->MajorGrid->IsVisible = true;
				 if (fabs(chart->Get_axis_scale_x()[0] - NA) > 1e-3)
					 myPane->XAxis->Scale->Min = chart->Get_axis_scale_x()[0];
				 else
					 myPane->XAxis->Scale->MinAuto = true;
				 if (fabs(chart->Get_axis_scale_x()[1] - NA) > 1e-3)
					 myPane->XAxis->Scale->Max = chart->Get_axis_scale_x()[1];
				 else
					 myPane->XAxis->Scale->MaxAuto = true;
				 if (fabs(chart->Get_axis_scale_x()[2] - NA) > 1e-3)
					 myPane->XAxis->Scale->MajorStep = chart->Get_axis_scale_x()[2];
				 else
					 myPane->XAxis->Scale->MajorStepAuto = true;
				 if (fabs(chart->Get_axis_scale_x()[3] - NA) > 1e-3)
				 {
					 myPane->XAxis->Scale->MinorStep = chart->Get_axis_scale_x()[3];
					 if (chart->Get_axis_scale_x()[3] == 0.0)
						 // remove minor tics
						 myPane->XAxis->MinorTic->Size = 0;
				 }
				 else
					 myPane->XAxis->Scale->MinorStepAuto = true;
				 if (chart->Get_axis_scale_x()[4] == 10.0)
					 myPane->XAxis->Type = AxisType::Log;

				 // Make the Y axis scale red
				 // myPane->YAxis->Scale->FontSpec->FontColor = Color::Red;
				 // myPane->YAxis->Title->FontSpec->FontColor = Color::Red;
				 // turn off the opposite tics so the Y tics don`t show up on the Y2 axis
				 if (Y2show)
				 {
					 myPane->YAxis->MajorTic->IsOpposite = false;
					 myPane->YAxis->MinorTic->IsOpposite = false;
				 }
				 // Don`t display the Y zero line
				 myPane->YAxis->MajorGrid->IsZeroLine = false;
				 // Align the Y axis labels so they are flush to the axis
				 myPane->YAxis->Scale->Align = AlignP::Inside;
				 myPane->YAxis->MajorGrid->IsVisible = true;
				 if (fabs(chart->Get_axis_scale_y()[0] - NA) > 1e-3)
					 myPane->YAxis->Scale->Min = chart->Get_axis_scale_y()[0];
				 else
					 myPane->YAxis->Scale->MinAuto = true;
				 if (fabs(chart->Get_axis_scale_y()[1] - NA) > 1e-3)
					 myPane->YAxis->Scale->Max = chart->Get_axis_scale_y()[1];
				 else
					 myPane->YAxis->Scale->MaxAuto = true;
				 if (fabs(chart->Get_axis_scale_y()[2] - NA) > 1e-3)
					 myPane->YAxis->Scale->MajorStep = chart->Get_axis_scale_y()[2];
				 else
					 myPane->YAxis->Scale->MajorStepAuto = true;
				 if (fabs(chart->Get_axis_scale_y()[3] - NA) > 1e-3)
				 {
					 myPane->YAxis->Scale->MinorStep = chart->Get_axis_scale_y()[3];
					 if (chart->Get_axis_scale_y()[3] == 0.0)
						 // remove minor tics
						 myPane->YAxis->MinorTic->Size = 0;
				 }
				 else
					 myPane->YAxis->Scale->MinorStepAuto = true;
				 if (chart->Get_axis_scale_y()[4] == 10.0)
					 myPane->YAxis->Type = AxisType::Log;

				 // Enable the Y2 axis display
				 if (Y2show)
				 {
					 myPane->Y2Axis->IsVisible = true;
					 // Make the Y2 axis scale blue
					 // myPane->Y2Axis->Scale->FontSpec->FontColor = Color::Blue;
					 // myPane->Y2Axis->Title->FontSpec->FontColor = Color::Blue;
					 // turn off the opposite tics so the Y2 tics don`t show up on the Y axis
					 myPane->Y2Axis->MajorTic->IsOpposite = false;
					 myPane->Y2Axis->MinorTic->IsOpposite = false;
					 // Don`t display the Y2 axis grid lines
					 myPane->Y2Axis->MajorGrid->IsVisible = false;
					 // Align the Y2 axis labels so they are flush to the axis
					 myPane->Y2Axis->Scale->Align = AlignP::Inside;

					 if (fabs(chart->Get_axis_scale_y2()[0] - NA) > 1e-3)
						 myPane->Y2Axis->Scale->Min = chart->Get_axis_scale_y2()[0];
					 else
						 myPane->Y2Axis->Scale->MinAuto = true;
					 if (fabs(chart->Get_axis_scale_y2()[1] - NA) > 1e-3)
						 myPane->Y2Axis->Scale->Max = chart->Get_axis_scale_y2()[1];
					 else
						 myPane->Y2Axis->Scale->MaxAuto = true;
					 if (fabs(chart->Get_axis_scale_y2()[2] - NA) > 1e-3)
						 myPane->Y2Axis->Scale->MajorStep = chart->Get_axis_scale_y2()[2];
					 else
						 myPane->Y2Axis->Scale->MajorStepAuto = true;
					 if (fabs(chart->Get_axis_scale_y2()[3] - NA) > 1e-3)
					 {
						 myPane->Y2Axis->Scale->MinorStep = chart->Get_axis_scale_y2()[3];
						 if (chart->Get_axis_scale_y2()[3] == 0.0)
							 // remove minor tics
							 myPane->Y2Axis->MinorTic->Size = 0;
					 }
					 else
						 myPane->Y2Axis->Scale->MinorStepAuto = true;
					 if (chart->Get_axis_scale_y2()[4] == 10.0)
						 myPane->Y2Axis->Type = AxisType::Log;
				 }

				 myPane->XAxis->MinorTic->IsOutside = false;
				 myPane->XAxis->MajorTic->IsOutside = false;
				 myPane->YAxis->MinorTic->IsOutside = false;
				 myPane->YAxis->MajorTic->IsOutside = false;
				 myPane->Y2Axis->MinorTic->IsOutside = false;
				 myPane->Y2Axis->MajorTic->IsOutside = false;

				 // Fill the axis background with a gradient
				 //myPane->Chart->Fill = gcnew Fill( Color::White, Color::LightYellow, 45.0f ); /* FromArgb(255, 255, 224) */
				if (this->background)
				{
					myPane->Chart->Fill = gcnew Fill( Color::White, Color::FromArgb(255, 255, 230), 45.0f );
				}
				else
				{
					myPane->Chart->Fill = gcnew Fill( Color::White, Color::White, 45.0f );
				}
				if (this->grid)
				{
					myPane->XAxis->MajorGrid->IsVisible = true;
					myPane->YAxis->MajorGrid->IsVisible = true;
				}
				else
				{
					myPane->XAxis->MajorGrid->IsVisible = false;
					myPane->YAxis->MajorGrid->IsVisible = false;
				}
				 // normalize pane size...
				 myPane->BaseDimension = 8.0F;
				 // increase bottom margin to accommodate text options...
				 myPane->Margin->Bottom = 15.0F;

				 // Make sure auto scale, Refresh
				 zg1->AxisChange();
				 zg1->Refresh();

			 }

	public: void CreateGraph( ZedGraphControl ^z1 )	{
				// Get a reference to the GraphPane instance in the ZedGraphControl
				GraphPane ^myPane = z1->GraphPane;

				// lock thread
				while (0 != System::Threading::Interlocked::CompareExchange(this->chartobject_ptr->usingResource, 2, 0))
				{
					System::Threading::Thread::Sleep(5);
				}

				try
				{

					DefineCurves(myPane, 0);

					// Add text boxes with instructions...
					GOL_no_hints = gcnew ZedGraph::GraphObjList(myPane->GraphObjList);

					TextObj ^text;
					text = gcnew TextObj(
						L" Click right mouse for options... \0",
						0.01f, 0.99f, CoordType::PaneFraction, AlignH::Left, AlignV::Bottom );
					text->FontSpec->StringAlignment = StringAlignment::Near;
					text->FontSpec->Size = 10;
					text->FontSpec->FontColor = Color::Red;
					text->ZOrder = ZOrder::H_BehindAll;
					myPane->GraphObjList->Add( text );
					text = gcnew TextObj(
						L" Press Alt + F4 to quit",
						0.81f, 0.99f, CoordType::PaneFraction, AlignH::Left, AlignV::Bottom );
					text->FontSpec->StringAlignment = StringAlignment::Near;
					text->FontSpec->Size = 10;
					text->FontSpec->FontColor = Color::Red;
					text->ZOrder = ZOrder::H_BehindAll;
					myPane->GraphObjList->Add( text );

					GOL_hints = gcnew ZedGraph::GraphObjList(myPane->GraphObjList);
					if (this->hints)
					{
						myPane->GraphObjList = GOL_hints;
					}
					else
					{
						myPane->GraphObjList = GOL_no_hints;
					}

					// Enable scrollbars if needed...
					/*z1->IsShowHScrollBar = true;
					z1->IsShowVScrollBar = true;
					z1->IsAutoScrollRange = true;
					z1->IsScrollY2 = true;*/

					// OPTIONAL: Show tooltips when the mouse hovers over a point
					z1->IsShowPointValues = false;
					z1->PointValueEvent += gcnew ZedGraphControl::PointValueHandler( this,
						&Form1::MyPointValueHandler );

					// OPTIONAL: Add a custom context menu item
					z1->ContextMenuBuilder += gcnew	ZedGraphControl::ContextMenuBuilderEventHandler(
						this, &Form1::MyContextMenuBuilder );

					// OPTIONAL: Handle the Zoom Event
					z1->ZoomEvent += gcnew ZedGraphControl::ZoomEventHandler( this,
						&Form1::MyZoomEvent );

					// Size the control to fit the window
					SetSize();

					// Tell ZedGraph to calculate the axis ranges
					// Note that you MUST call this after enabling IsAutoScrollRange, since AxisChange() sets
					// up the proper scrolling parameters

					z1->AxisChange();
					// Make sure the Graph gets redrawn
					z1->Invalidate();
					timer1->Interval = this->chartobject_ptr->Get_update_time_chart();
					timer1->Enabled = true;
					timer1->Start();

					tickStart = Environment::TickCount;
				}
				catch (...)
				{
					//unlock thread
					int n = System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
					assert(n == 2);
					throw;
				}

				//unlock thread
				int n = System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
				assert(n == 2);
			}

			/// <summary>
			/// Display customized tooltips when the mouse hovers over a point
			/// </summary>
			System::String ^MyPointValueHandler( ZedGraphControl ^control, GraphPane ^pane,
				CurveItem ^curve, int iPt ) {
					// Get the PointPair that is under the mouse
					PointPair pt = curve[iPt];
					return curve->Label->Text + " is " + pt.Y.ToString( "e3" ) + " units at X = " + pt.X.ToString( "e3" );
			}

			// Add some explanation to the menu..
			void MyContextMenuBuilder( ZedGraphControl ^control,
				System::Windows::Forms::ContextMenuStrip ^menuStrip,
				Point mousePt,
				ZedGraphControl::ContextMenuObjectState objState ) {

					// removes Copy
					menuStrip->Items->RemoveAt(0);
					// removes Save Image As
					menuStrip->Items->RemoveAt(0);

					ToolStripMenuItem ^item = gcnew ToolStripMenuItem();
					item->Text = L"Zoom: left mouse + drag\nPan: middle mouse + drag";
					menuStrip->Items->Insert(5, item );

					ToolStripMenuItem ^item3 = gcnew ToolStripMenuItem();
					item3->Text = L"Chart options...";
					item3->Click += gcnew System::EventHandler(this, &zdg_ui2::Form1::SetChartOptions );
					menuStrip->Items->Insert(0, item3 );

					//ToolStripMenuItem ^item5 = gcnew ToolStripMenuItem();
					//item5->Text = L"Toggle Hints";
					//item5->Click += gcnew System::EventHandler(this, &zdg_ui2::Form1::ToggleHints );
					//menuStrip->Items->Insert(0, item5 );

					ToolStripMenuItem ^item2 = gcnew ToolStripMenuItem();
					item2->Text = L"Save Data to File...";
					item2->Click += gcnew System::EventHandler(this, &zdg_ui2::Form1::SaveCurves );
					menuStrip->Items->Insert(0, item2 );

					ToolStripMenuItem ^item4 = gcnew ToolStripMenuItem();
					item4->Text = L"Save Image As...";
					item4->Click += gcnew System::EventHandler(this, &zdg_ui2::Form1::SaveImage );
					menuStrip->Items->Insert(0, item4 );

			}

			void form_error_msg( std::string estring )
			{
				if (this->chartobject_ptr != NULL)
				{
					if (Phreeqc* ptr = this->chartobject_ptr->Get_phreeqc())
					{
						ptr->error_msg(estring.c_str(), CONTINUE);
					}
				}
				else
				{
					std::cerr << "ERROR: " << estring << std::endl;
				}
			}

			void SaveCurves( System::Object ^sender, System::EventArgs ^e )
			{
				SaveFileDialog^ saveFileDialog1 = gcnew SaveFileDialog;
				//TCHAR dir[MAX_PATH];
				//::GetCurrentDirectory(MAX_PATH, dir);
				//String ^d = gcnew String(dir);
				//saveFileDialog1->InitialDirectory = d;
				saveFileDialog1->FileName = "curves.u_g";
				saveFileDialog1->Filter = "User graph files (*.u_g)|*.u_g|txt files (*.txt)|*.txt|All files (*.*)|*.*";
				saveFileDialog1->FilterIndex = 1;
				saveFileDialog1->RestoreDirectory = true;
#undef OK 
				if ( saveFileDialog1->ShowDialog() == System::Windows::Forms::DialogResult::OK )
				{
					std::string file_name;
					ToString(saveFileDialog1->FileName, file_name);
					std::ofstream f_out(file_name.c_str(), std::ifstream::out);

					if (!f_out.is_open())
					{
						std::ostringstream estream;
						estream << "Could not open file to save curves for USER_GRAPH " << file_name;
						form_error_msg(estream.str());
						return;
					}

					// write headings
					size_t max_points = 0; 
					f_out.precision(4);
					for (int i = 0; i < zg1->GraphPane->CurveList->Count; i++) 
					{
						LineItem  ^curve;
						curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
						// Get the PointPairList
						IPointListEdit  ^ip = (IPointListEdit^) curve->Points;

						// Calculate max_points
						if ((size_t) ip->Count > max_points)
							max_points = ip->Count;

						// write headers
						std::string s_std;
						ToString(curve->Label->Text, s_std);
						f_out.width(12);
						f_out << "x" << "\t";
						f_out.width(12);
						if (s_std.size() > 0) 
						{
							f_out << s_std << "\t";
						}
						else
						{
							f_out << "y" << "\t";
						}
					}

					f_out << std::endl;

					// write data
					size_t i2 = 0;
					f_out << std::scientific;
					f_out.precision(4);

					while (i2 < max_points)
					{
						for (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
						{
							LineItem  ^curve;
							curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
							// Get the PointPairList
							IPointListEdit  ^ip = (IPointListEdit^) curve->Points;					
							if (i2 < (size_t) ip->Count)
							{
								//double x = ip[i]->X;
								f_out.width(12);
								f_out << ip[(int) i2]->X << "\t";
								f_out.width(12);
								f_out << ip[(int) i2]->Y << "\t";
							}
							else if (i2 < max_points)
							{
								f_out.width(13);
								f_out << "\t";
								f_out.width(13);
								f_out << "\t";
							}
						}
						f_out << std::endl;
						i2++;
					}
					f_out.close();


				}
				else
				{
					// no dialog
					std::ostringstream estream;
					estream << "Could not open dialog to save curves. ";
					form_error_msg(estream.str());
					return;
				}
				return;
			}

			// Respond to a Zoom Event
			void MyZoomEvent( ZedGraphControl ^control, ZoomState ^oldState, ZoomState ^newState )
			{
				// Here we get notification everytime the user zooms
			}
			void SetChartOptions( System::Object ^sender, System::EventArgs ^e )
			{
				// Create form
				Form ^graphOptions = gcnew Form;
				graphOptions->Text = "Chart options";
				graphOptions->ShowIcon = false;
				graphOptions->BringToFront();
				graphOptions->MinimumSize = System::Drawing::Size(255, 230);
				graphOptions->MaximumSize = System::Drawing::Size(255, 230);
				graphOptions->Size = System::Drawing::Size(255, 230);
				graphOptions->StartPosition = System::Windows::Forms::FormStartPosition::CenterParent;

				// Check box for hints
				CheckBox ^cb1 = gcnew CheckBox;
				cb1->Appearance = Appearance::Normal;
				cb1->ThreeState = false;
				cb1->AutoCheck = true;
				cb1->Location = System::Drawing::Point(5, 10);
				if (this->hints)
				{
					cb1->CheckState = CheckState::Checked;
				}
				else
				{
					cb1->CheckState = CheckState::Unchecked;
				}
				cb1->Text = "Show hints";
				cb1->Visible = true;
				graphOptions->Controls->Add(cb1);

				// Check box for background color
				CheckBox ^cb2 = gcnew CheckBox;
				cb2->Appearance = Appearance::Normal;
				cb2->ThreeState = false;
				cb2->AutoCheck = true;
				cb2->Location = System::Drawing::Point(5, 60);

				if (this->background)
				{
					cb2->CheckState = CheckState::Checked;
				}
				else
				{
					cb2->CheckState = CheckState::Unchecked;
				}
				cb2->Text = "Show colored background";
				cb2->AutoSize = true;
				cb2->Visible = true;
				graphOptions->Controls->Add(cb2);

				// checkbox for grid
				CheckBox ^cb3 = gcnew CheckBox;
				cb3->Appearance = Appearance::Normal;
				cb3->ThreeState = false;
				cb3->AutoCheck = true;
				cb3->Location = System::Drawing::Point(5, 110);

				if (this->grid)
				{
					cb3->CheckState = CheckState::Checked;
				}
				else
				{
					cb3->CheckState = CheckState::Unchecked;
				}
				cb3->Text = "Show grid lines";
				cb3->Visible = true;
				graphOptions->Controls->Add(cb3);

				// done button for Form
				Button^ button1 = gcnew Button;
				button1->DialogResult = System::Windows::Forms::DialogResult::OK;
				button1->Text = "Done";
				button1->Location = System::Drawing::Point(75, 160);
				graphOptions->Controls->Add(button1);
				graphOptions->AcceptButton = button1;

				// cancel button for Form
				Button^ button2 = gcnew Button;
				button2->DialogResult = System::Windows::Forms::DialogResult::Cancel;
				button2->Text = "Cancel";
				button2->Location = System::Drawing::Point(155, 160);
				graphOptions->Controls->Add(button2);
				graphOptions->CancelButton = button2;


				if (graphOptions->ShowDialog() == System::Windows::Forms::DialogResult::OK)
				{
					this->hints = (cb1->CheckState == CheckState::Checked);
					this->background = (cb2->CheckState == CheckState::Checked);
					this->grid = (cb3->CheckState == CheckState::Checked);
				}

				if (this->background)
				{
					zg1->GraphPane->Chart->Fill = gcnew Fill( Color::White, Color::FromArgb(255, 255, 230), 45.0f );
				}
				else
				{
					zg1->GraphPane->Chart->Fill = gcnew Fill( Color::White, Color::White, 45.0f );
				}
				if (this->grid)
				{
					zg1->GraphPane->XAxis->MajorGrid->IsVisible = true;
					zg1->GraphPane->YAxis->MajorGrid->IsVisible = true;
				}
				else
				{
					zg1->GraphPane->XAxis->MajorGrid->IsVisible = false;
					zg1->GraphPane->YAxis->MajorGrid->IsVisible = false;
				}
				if (this->hints)
				{
					zg1->GraphPane->GraphObjList = GOL_hints;
				}
				else
				{
					zg1->GraphPane->GraphObjList = GOL_no_hints;
				}
				zg1->Refresh();
			}
			void SaveImage( System::Object ^sender, System::EventArgs ^e )
			{
				ZedGraph::GraphObjList ^copy = gcnew ZedGraph::GraphObjList(zg1->GraphPane->GraphObjList);
				zg1->GraphPane->GraphObjList = GOL_no_hints;
				zg1->SaveAs();
				zg1->GraphPane->GraphObjList = copy;
			}
			void BatchSaveImage( )
			{
				ChartObject *chart = this->chartobject_ptr;
				assert(chart->Get_batch() > 0);
				// Save GraphObjList
				ZedGraph::GraphObjList ^GOL_copy = gcnew ZedGraph::GraphObjList(zg1->GraphPane->GraphObjList);

				// Don`t write red hint boxes
				zg1->GraphPane->GraphObjList = GOL_no_hints;

				// Set background
				if (chart->Get_batch_background())
				{
					zg1->GraphPane->Chart->Fill = gcnew Fill( Color::White, Color::FromArgb(255, 255, 230), 45.0f );
				}
				else
				{
					zg1->GraphPane->Chart->Fill = gcnew Fill( Color::White, Color::White, 45.0f );
				}
				// Set grid
				if (chart->Get_batch_grid())
				{
					zg1->GraphPane->XAxis->MajorGrid->IsVisible = true;
					zg1->GraphPane->YAxis->MajorGrid->IsVisible = true;
				}
				else
				{
					zg1->GraphPane->XAxis->MajorGrid->IsVisible = false;
					zg1->GraphPane->YAxis->MajorGrid->IsVisible = false;
				}
				// Save the graph
				if (this->zg1)
				{
					ImageFormat ^format = ImageFormat::Png;
					switch (chart->Get_batch())
					{
					case 2:
						format = ImageFormat::Png;
						break;
					case 3:
						format = ImageFormat::Gif;
						break;
					case 4:
						format = ImageFormat::Jpeg;
						break;
					case 5:
						format = ImageFormat::Tiff;
						break;
					case 6:
						format = ImageFormat::Bmp;
						break;
					default:
						break;
					}
					switch (chart->Get_batch())
					{
					case 1: // emf
						{
							System::String ^fn = gcnew System::String(chart->Get_batch_fn().c_str());
							System::Drawing::Imaging::Metafile ^metaFile = this->zg1->MasterPane->GetMetafile();
							metaFile->Save(fn);
						}
						break;
					case 2: // bitmaps
					case 3:
					case 4:
					case 5:
					case 6:
						{
							System::String ^fn = gcnew System::String(chart->Get_batch_fn().c_str());
							System::IO::FileStream ^myStream = gcnew System::IO::FileStream(fn, System::IO::FileMode::Create);
							zg1->MasterPane->GetImage()->Save( myStream, format);
							myStream->Close();
						}
					default:
						break;
					}
				}

			}
   private: void timer1_Tick(System::Object ^sender, System::EventArgs ^e )
			{
				LineItem  ^curve;
				ChartObject *chart = this->chartobject_ptr;
				if (chart == NULL) return;

				//lock for thread
				while (0 != System::Threading::Interlocked::CompareExchange(chart->usingResource, 3, 0))
				{
					System::Threading::Thread::Sleep(5);
				}

				try
				{
					if (this->chartobject_ptr->Get_curve_added())
					{
						DefineCurves(zg1->GraphPane, zg1->GraphPane->CurveList->Count);
					}
					else if (this->chartobject_ptr->Get_point_added())
					{

						// Make list of curves
						std::vector<CurveObject *> Curves; 
						size_t j;
						for (j = 0; j < chart->Get_CurvesCSV().size(); j++)
						{
							Curves.push_back(chart->Get_CurvesCSV()[j]);
						}
						for (j = 0; j < chart->Get_Curves().size(); j++)
						{
							Curves.push_back(chart->Get_Curves()[j]);
						}
						// Add points to curves ...
						for (int i = 0; i < zg1->GraphPane->CurveList->Count; i++) 
						{
							curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
							// Get the PointPairList
							IPointListEdit  ^ip = (IPointListEdit^) curve->Points;
							if ((size_t) ip->Count < Curves[i]->Get_x().size())
							{
								if (Curves[i]->Get_y_axis() == 2)
									Y2 = true;
								else
									Y2 = false;
								for ( size_t i2 = ip->Count; i2 < Curves[i]->Get_x().size(); i2++ )
								{
									if ((LogX && Curves[i]->Get_x()[i2] <=0)
										|| (LogY && !Y2 && Curves[i]->Get_y()[i2] <=0)
										|| (LogY2 && Y2 && Curves[i]->Get_y()[i2] <=0))
										continue;
									else
										ip->Add(Curves[i]->Get_x()[i2], Curves[i]->Get_y()[i2] );
								}
							}
						}
						// Add points to curves ...

						//size_t i, k;
						//k = 0;
						//for (i = 0; i < Curves.size(); i++) 
						//{
						//	if (Curves[i]->Get_x().size() == 0) continue;
						//	curve =  (LineItem ^) zg1->GraphPane->CurveList[k++];
						//	// Get the PointPairList
						//	IPointListEdit  ^ip = (IPointListEdit^) curve->Points;
						//	if ((size_t) ip->Count < Curves[i]->Get_x().size())
						//	{
						//		if (Curves[i]->Get_y_axis() == 2)
						//			Y2 = true;
						//		else
						//			Y2 = false;
						//		for ( size_t i2 = ip->Count; i2 < Curves[i]->Get_x().size(); i2++ )
						//		{
						//			if ((LogX && Curves[i]->Get_x()[i2] <=0)
						//				|| (LogY && !Y2 && Curves[i]->Get_y()[i2] <=0)
						//				|| (LogY2 && Y2 && Curves[i]->Get_y()[i2] <=0))
						//				continue;
						//			else
						//				ip->Add(Curves[i]->Get_x()[i2], Curves[i]->Get_y()[i2] );
						//		}
						//	}
						//}
						/* explicitly reset the max in case of log scale, zedgraphs doesn`t do this... */
						if ((fabs(chart->Get_axis_scale_x()[1] - NA) < 1e-3) && zg1->GraphPane->XAxis->Type == AxisType::Log)
						{
							double max = -1e99;
							for  (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
							{
								if (Curves[i]->Get_x()[Curves[i]->Get_x().size() - 1] > max)
									max = Curves[i]->Get_x()[Curves[i]->Get_x().size() - 1];
							}
							max += pow(10.0, log10(max / 3));
							zg1->GraphPane->XAxis->Scale->Max = max;
						}
						if ((fabs(chart->Get_axis_scale_y()[1] - NA) < 1e-3) && zg1->GraphPane->YAxis->Type == AxisType::Log)
						{
							double max = -1e99;
							for  (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
							{
								curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
								if (curve->IsY2Axis) continue;
								if (Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1] > max)
									max = Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1];
							}
							max += pow(10.0, log10(max / 3));
							zg1->GraphPane->YAxis->Scale->Max = max;
						}
						if ((fabs(chart->Get_axis_scale_y2()[1] - NA) < 1e-3) && zg1->GraphPane->Y2Axis->Type == AxisType::Log)
						{
							double max = -1e99;
							for  (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
							{
								curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
								if (!curve->IsY2Axis) continue;
								if (Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1] > max)
									max = Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1];
							}
							max += pow(10.0, log10(max / 3));
							zg1->GraphPane->Y2Axis->Scale->Max = max;
						}
						zg1->AxisChange();
						zg1->Refresh();
					}
					//
					// Following asserts may not be true for Log scales 
					// negative values are rejected, so chart may have fewer points than phreeqc
					//
					//for (size_t j = 0; j < chart->Get_CurvesCSV().size(); j++)
					//{
					//	if (zg1->GraphPane->CurveList[j]->Points->Count != chart->Get_CurvesCSV()[j]->Get_x().size())
					//	{
					//		fprintf(stderr, "graph points = %d\n", zg1->GraphPane->CurveList[j]->Points->Count);
					//		fprintf(stderr, "phreeqc points = %d\n", chart->Get_CurvesCSV()[j]->Get_x().size());
					//	}
					//	assert(zg1->GraphPane->CurveList[j]->Points->Count == chart->Get_CurvesCSV()[j]->Get_x().size());
					//}
					//for (int j = chart->Get_CurvesCSV().size(); j < zg1->GraphPane->CurveList->Count; j++) 
					//{
					//	int k = j - chart->Get_CurvesCSV().size();
					//	if (zg1->GraphPane->CurveList[j]->Points->Count != chart->Get_Curves()[k]->Get_x().size())
					//	{
					//		fprintf(stderr, "%d %d graph points = %d\n", j, k, zg1->GraphPane->CurveList[j]->Points->Count);
					//		fprintf(stderr, "phreeqc points = %d\n", chart->Get_Curves()[k]->Get_x().size());
					//	}
					//	assert(zg1->GraphPane->CurveList[j]->Points->Count == chart->Get_Curves()[k]->Get_x().size());
					//}
					chart->Set_point_added(false);
					if (chart->Get_end_timer())
					{
						timer1->Stop();
						chart->Set_done(true);
						phreeqc_done = true;

						int batch = chart->Get_batch();
						if (batch > 0)
						{
							BatchSaveImage();
						}

						//unlock thread before setting chartobject_ptr to NULL
						int n = System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
						assert(n == 3);

						//this->phreeqc_ptr = NULL;
						this->chartobject_ptr = NULL;
						if (batch >= 0)
						{
							this->Close();
						}
						return;
					}
				}
				catch(...)
				{
					int n = System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
					assert(n == 3);

					throw;
				}

				//unlock thread
				int n = System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
				assert(n == 3);
				//tickStart = Environment::TickCount;
				return;
			}

			 void ToString(System::String^ src, std::string& dest)
			 {
				 using namespace System::Runtime::InteropServices;
				 const char* chars = (const char*)(Marshal::StringToHGlobalAnsi(src)).ToPointer();
				 dest = chars;
				 Marshal::FreeHGlobal(IntPtr((void*)chars));
			 }
			 ~Form1() {
				 if (this->zg1) delete zg1;
				 //if (this->timer1) delete timer1);
				 if (components) {
					 delete components;
				 }
			 }
	public: ZedGraph::ZedGraphControl ^zg1;
	private: System::Windows::Forms::Timer ^timer1;
	private: System::ComponentModel::Container ^components;
			 ChartObject * chartobject_ptr;

	public:
#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent()
		{
			//System::Threading::Thread::Sleep(5000);
			this->components = (gcnew System::ComponentModel::Container());
			this->zg1 = (gcnew ZedGraph::ZedGraphControl());
			this->timer1 = (gcnew System::Windows::Forms::Timer( this->components ));
			this->SuspendLayout();
			// 
			// zg1
			// 
			this->zg1->Location = System::Drawing::Point(12, 12);
			this->zg1->Name = L"zg1";
			this->zg1->ScrollGrace = 0;
			this->zg1->ScrollMaxX = 0;
			this->zg1->ScrollMaxY = 0;
			this->zg1->ScrollMaxY2 = 0;
			this->zg1->ScrollMinX = 0;
			this->zg1->ScrollMinY = 0;
			this->zg1->ScrollMinY2 = 0;
			this->zg1->Size = System::Drawing::Size(this->chartobject_ptr->Get_PanelWidth() - 2 * 12, chartobject_ptr->Get_PanelHeight() - 2 * 12);
			this->zg1->TabIndex = 0;
			this->timer1->Tick += gcnew System::EventHandler( this, &Form1::timer1_Tick );
			// 
			// Form1
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->AutoValidate = System::Windows::Forms::AutoValidate::EnablePreventFocusChange;
			this->ClientSize = System::Drawing::Size(this->chartobject_ptr->Get_PanelWidth(), chartobject_ptr->Get_PanelHeight());
			this->Controls->Add(this->zg1);
			this->Name = L"Form1";
			this->StartPosition = System::Windows::Forms::FormStartPosition::WindowsDefaultLocation;//:CenterScreen;
			//this->Text = L"PHREEQC chart";
			System::String ^desc = gcnew String(this->chartobject_ptr->Get_description().c_str());
			this->Text = L"Phreeqc USER_GRAPH " + this->chartobject_ptr->Get_n_user() + " " + desc;
			this->TopMost = false;
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(Form1::typeid));
			try
			{
				//this->Icon = gcnew System::Drawing::Icon("c:\\phreeqc\\phreex.ico");
				this->Icon = (cli::safe_cast<System::Drawing::Icon^>(resources->GetObject(L"$this.Icon")));
			}
			catch (...)
			{
			}

			this->FormClosing += gcnew System::Windows::Forms::FormClosingEventHandler(this, &Form1::MyFormClosingEventHandler);
			this->Load += gcnew System::EventHandler(this, &Form1::Form1_Load);
			this->Resize += gcnew System::EventHandler(this, &Form1::Form1_Resize);
			this->ResumeLayout(false);
		}
#pragma endregion
	};
}