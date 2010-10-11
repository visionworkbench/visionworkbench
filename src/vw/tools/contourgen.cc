// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <list>
#include <map>
#include <algorithm>
#include <utility>

#include <cairomm/context.h>
#include <cairomm/surface.h>

//#include <boost/algorithm/minmax_element.hpp>
#include <boost/program_options.hpp>

#include <vw/Core/Log.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Image.h>
#include <vw/FileIO.h>

#include <cmath>

#include "contour.h"

namespace po = boost::program_options;

// fwd declaration of function in FitCurves.cpp
BezierContour FitCurve(PointContour &contour, double error);

vw::ImageView<float> load_dem_from_file(std::string file_in) {
    //    Open the file for reading
    vw::ImageView<float> dem;
    vw::DiskImageResource *r = vw::DiskImageResource::open( file_in );
    vw::vw_out(vw::InfoMessage) << "Loading DEM from file" << std::endl;
    vw::vw_out(vw::DebugMessage) << r->rows()
        << "x" << r->cols() << "x" << r->planes()
        << "  " << r->channels() << " channel(s)"
        << " type: " << r->type()
        << " pxfmt: " << r->pixel_format()
        << " chtype: " << r->channel_type()
        << "\n";
    //    Read the data
    //tpc.set_progress_text("Status (loading image):   " );
    vw::read_image(dem, *r);
    delete r;
    return dem;
}

void read_segments_from_file(SegmentList &segment_list,
        std::string file_in, int &rows, int &cols)
{
    vw::vw_out(vw::InfoMessage)
        << "Reading contour segments from file: " << file_in << "\t" << std::endl;
    std::ifstream ifs;
    ifs.open(file_in.c_str(), std::ifstream::in);

    // first line should have image dimensions
    if (ifs.good()) ifs >> std::skipws >> rows >> cols;

    // rest of file should have segment coordinates and level
    while (ifs.good()) {
        float ax,ay,bx,by;
        int level;
        ContourSegment seg;
        ifs >> std::skipws >> ax >> ay >> bx >> by >> level;
        seg.a = ContourPoint(ax,ay);
        seg.b = ContourPoint(bx,by);
        seg.level = level;
        segment_list.push_back(seg);
    }
    ifs.close();
    vw::vw_out(vw::DebugMessage)
        << "\tRead " << segment_list.size() << " segments" << std::endl;
}

void load_segments_into_pcs(SegmentList &segment_list, PointContourSet &cset) {
    SegmentList::iterator seglist_iter;

    // Load segments into PointContourSet
    vw::vw_out(vw::InfoMessage)
        << "Loading segments into PointContourSet" << std::endl;

    for (seglist_iter = segment_list.begin();
         seglist_iter != segment_list.end();
         seglist_iter++) {
        add_segment(cset, *seglist_iter);
    }
}

void fit_Bezier_curves_to_contours(PointContourSet &cset,
        BezierContourSet &bcset, float error)
{
    vw::vw_out(vw::InfoMessage)
        << "Fitting Bezier curves to contours" << std::endl;
    PointContourSet::iterator cset_iter;
    PointContour c;
    int level;

    for (cset_iter = cset.begin(); cset_iter != cset.end(); cset_iter++) {
        level = (*cset_iter).first;
        c = (*cset_iter).second;
        BezierContour bc = FitCurve(c, error);
        bcset.insert(make_pair(level, bc));
        /*
        std::cout << "adding contour: " << bc.size() << " elems" << std::endl;
        for (BezierContour::iterator bciter = bc.begin(); bciter != bc.end(); ++bciter) {
            for (int i = 0; i < 4; i++) {
                std::cout << "\t" << (*bciter)[i] << std::endl;
            }
        }
        */
    }
}

Cairo::RefPtr<Cairo::Context> create_png_surface(int dem_w, int dem_h,
        std::string image_file)
{
    vw::vw_out(vw::InfoMessage) << "Setting up PNG output surface" << std::endl;

    Cairo::RefPtr<Cairo::ImageSurface> surface =
        Cairo::ImageSurface::create_from_png(image_file);
    Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);

    int img_w = surface->get_width();
    int img_h = surface->get_height();
    double scale_w = double(img_w)/double(dem_w);
    double scale_h = double(img_h)/double(dem_h);
    cr->scale(scale_w, scale_h);

    vw::vw_out(vw::DebugMessage) <<
        "\tDEM dimensions: " << "(" << dem_w << "x" << dem_h << ")" << std::endl
        << "\tImg dimensions: " << "(" << img_w << "x" << img_h << ")" << std::endl
        << "\tScale factors:  " << "(" << scale_w << "x" << scale_h << ")" << std::endl;

    cr->save();
    return cr;
}

Cairo::RefPtr<Cairo::Context> create_svg_surface(int width_out, int height_out,
        std::string file_out)
{
    vw::vw_out(vw::InfoMessage) << "Setting up SVG output surface " << std::endl;

    Cairo::RefPtr<Cairo::SvgSurface> surface =
        Cairo::SvgSurface::create(file_out, width_out, height_out);
    Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);

    cr->save();
    cr->set_source_rgb(1.0, 1.0, 1.0);
    cr->paint();
    //cr->restore();
    cr->save();
    return cr;
}

Cairo::RefPtr<Cairo::Context>
create_output_surface(std::string output_type, int cols, int rows,
        std::string image_file, std::string file_out)
{
    Cairo::RefPtr<Cairo::Context> cr;
    if (output_type == "png") {
        // Create PNG surface
        cr = create_png_surface(cols, rows, image_file);
    } else {
        // Create SVG surface
        cr = create_svg_surface(cols, rows, file_out);
    }
    return cr;
}

void write_output_file(Cairo::RefPtr<Cairo::Context> cr,
        std::string output_type, std::string file_out)
{
    vw::vw_out(vw::InfoMessage) << "Writing output file: " << file_out << std::endl;

    if (output_type == "png") {
        cr->get_target()->write_to_png(file_out);
    } else {
        cr->show_page();
    }
}

void draw_point_contours(PointContourSet cset,
        Cairo::RefPtr<Cairo::Context> cr, float nodataval)
{
    vw::vw_out(vw::InfoMessage) << "Writing point contours to output surface\n";
    PointContourSet::iterator cset_iter;
    double line_width, tmp;
    line_width = 2.0; tmp = 0.0;
    cr->device_to_user_distance(line_width, tmp);
    cr->set_line_width(line_width);
    cr->set_source_rgb(1.0, 0.0, 0.0);

    int level = int(nodataval);
    for (cset_iter = cset.begin(); cset_iter != cset.end(); cset_iter++) {

        int newlevel = (*cset_iter).first;
        if (newlevel != level) {
            level = newlevel;
        }

        PointContour c = (*cset_iter).second;
        PointContour::iterator c_iter = c.begin();
        cr->begin_new_path();
        cr->move_to((*c_iter)[0],(*c_iter)[1]);
        while (++c_iter != c.end()) {
            cr->line_to((*c_iter)[0], (*c_iter)[1]);
        }
        cr->stroke();
    }
}

void draw_Bezier_contours(BezierContourSet bcset, Cairo::RefPtr<Cairo::Context> cr, float nodataval) {
    vw::vw_out(vw::InfoMessage) << "Writing Bezier contours to output surface\n";
    BezierContourSet::iterator bcset_iter;
    double line_width, tmp;
    line_width = 2.0; tmp = 0.0;
    //cr->device_to_user_distance(line_width, tmp);
    cr->set_line_width(line_width);
    cr->set_source_rgb(0.0, 0.0, 1.0);

    int level = int(nodataval);
    for (bcset_iter = bcset.begin(); bcset_iter != bcset.end(); bcset_iter++) {
        int newlevel = (*bcset_iter).first;
        if (newlevel != level) {
            level = newlevel;
        }

        BezierContour contour = (*bcset_iter).second;
        BezierContour::iterator c_iter = contour.begin();
        cr->begin_new_path();
        cr->move_to((*c_iter)[0][0], (*c_iter)[0][1]);
        while (++c_iter != contour.end()) {
            cr->curve_to((*c_iter)[1][0], (*c_iter)[1][1],
                         (*c_iter)[2][0], (*c_iter)[2][1],
                         (*c_iter)[3][0], (*c_iter)[3][1]);
        }
        cr->stroke();
    }
}



void write_points_to_file(std::string file_out, SegmentList segment_list, int rows, int cols) {
    vw::vw_out(vw::InfoMessage) << "Writing contour points to text file\n";
    SegmentList::iterator iter;
    std::ofstream ofs;
    ofs.open(file_out.c_str());
    ofs << rows << "\t" << cols << std::endl;

    for (iter = segment_list.begin(); iter != segment_list.end(); iter++) {
        ContourSegment s = *iter;
        ofs << std::setprecision(8) << std::fixed
            << s.a[0] << "\t" << s.a[1] << "\t" << s.b[0] << "\t"
            << s.b[1] << "\t" << s.level << std::endl;
    }
    ofs.close();
}


int main(int argc, char *argv[])
{
    bool verbose = false;
    std::string file_in;
    std::string file_out;

    // Parse command line options
    po::options_description opts("Allowed options");
    opts.add_options()
        ("help,?", "show help message")
        ("verbose,v", po::bool_switch(&verbose), "set verbose output")
        ("input-type,i", po::value<std::string>()->default_value("tiff"),
            "choose input type")
        ("image-file,m", po::value<std::string>()->default_value(""),
            "image file to write contours on")
        ("output-type,o", po::value<std::string>()->default_value("svg"),
            "choose output type")
        ("no-data-value,n", po::value<float>()->default_value(-10000.0),
            "set \"NO DATA\" value")
        ("contour-interval,c", po::value<int>()->default_value(100),
            "set contour interval")
    ;

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-file", po::value<std::string>(), "input file")
        ("output-file", po::value<std::string>(), "output file")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(opts).add(hidden);

    po::positional_options_description p;
    p.add("input-file", 1);
    p.add("output-file", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(cmdline_options).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << opts << std::endl;
        return 1;
    }

    if (vm.count("input-file")) {
        file_in = vm["input-file"].as<std::string>();
    } else {
        vw::vw_out(vw::ErrorMessage)
            << "ERROR: Must specify an input file" << std::endl;
        return 1;
    }

    if (vm.count("output-file")) {
        file_out = vm["output-file"].as<std::string>();
    } else {
        vw::vw_out(vw::ErrorMessage)
            << "ERROR: Must specify an output file" << std::endl;
        return 1;
    }

    std::string input_type = vm["input-type"].as<std::string>();
    std::string image_file = vm["image-file"].as<std::string>();
    std::string output_type = vm["output-type"].as<std::string>();
    float nodataval = vm["no-data-value"].as<float>();
    int cint = vm["contour-interval"].as<int>();

    vw::vw_log().console_log().rule_set().clear();
    vw::vw_log().console_log().rule_set().add_rule(vw::WarningMessage, "console");
    if (verbose) {
        vw::vw_log().console_log().rule_set().add_rule(vw::DebugMessage, "console");
    }

    if (output_type != "text" && output_type != "svg" && output_type != "png") {
        vw::vw_out(vw::ErrorMessage)
            << "ERROR: Output types other than \"text\", \"png\" or \"svg\" "
            << "not currently supported." << std::endl;
        return 1;
    }

    if (input_type != "tiff" && input_type != "text") {
        vw::vw_out(vw::ErrorMessage)
            << "ERROR: Input types other than \"tiff\" or \"text\" "
            << "not currently supported." << std::endl;
        return 1;
    }

    vw::vw_out(vw::DebugMessage) << "Input file: " << file_in << std::endl;
    vw::vw_out(vw::DebugMessage) << "Output file: " << file_out << std::endl;
    vw::vw_out(vw::DebugMessage) << "Input type: " << input_type << std::endl;
    vw::vw_out(vw::DebugMessage) << "Output type: " << output_type << std::endl;
    vw::vw_out(vw::DebugMessage) << "Image file: " << image_file << std::endl;
    vw::vw_out(vw::DebugMessage) << "No-data value: " << nodataval << std::endl;
    vw::vw_out(vw::DebugMessage) << "Contour interval: " << cint << std::endl;

    // done parsing options

    vw::ImageView<float> dem;
    //TerminalProgressCallback tpc();
    PointContourSet cset;
    SegmentList segment_list;
    int rows, cols;
    float error = 1.0e-3;

    if (input_type == "tiff") {
        // 1. Load DEM from file
        dem = load_dem_from_file(file_in);
        rows = dem.rows();
        cols = dem.cols();

        // 2. Run contouring algorithm on DEM
        conrec(dem, cset, cint, nodataval, segment_list);

    }
    else if (input_type == "text") {
        // read contour segments from file
        read_segments_from_file(segment_list, file_in, rows, cols);
    }

    if (output_type == "text") {
        write_points_to_file(file_out, segment_list, rows, cols);

    } else { // SVG or PNG
        // Load segments into point contour set
        load_segments_into_pcs(segment_list, cset);

        // Fit Bezier curves to contours
        BezierContourSet bcset;
        fit_Bezier_curves_to_contours(cset, bcset, error);

        Cairo::RefPtr<Cairo::Context> cr;

        cr = create_output_surface(output_type, cols, rows, image_file, file_out);

        //Draw contours onto surface
        //draw_point_contours(cset, cr, nodataval);
        draw_Bezier_contours(bcset, cr, nodataval);

        write_output_file(cr, output_type, file_out);

    }

    return 0;
}



