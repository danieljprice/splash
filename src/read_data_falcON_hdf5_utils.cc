// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file  read_falcON_utils.cc
///
/// \brief towards reading falcON HDF5 snapshots into splash
/// \date  12-Aug-2015
//
//  Copyright (C) 2015  Walter Dehnen.
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// PART 1:  interfaces with C-linkage
//          this part may be used like a header file
//
////////////////////////////////////////////////////////////////////////////////
extern "C" {
  /// maximum number of particles types expected
  /// (actually falcON currently only supports 3 particle types)
  static constexpr int max_num_types = 6;
  
  //
  //  1.1  routines to be provided by SPLASH and called from PART 2 below
  //

  ///
  /// transfer data to splash
  ///
  /// \param[in] icol   column of data
  /// \param[in] ndat   number of data == number of particles of type
  /// \param[in] data   array of ndat data
  /// \param[in] type   index of particle type
  ///
  /// \note called by read_falcON_snapshot()
  void read_falcON_data_into_splash(const int*icol, const int*ndat,
				    const double*data, const int*type);

  ///
  /// tell splash about the label of a data column
  ///
  /// \param[in] icol  column of data
  /// \param[in] name  label for this column
  ///
  /// \note called by read_falcON_snapshot()
  void set_splash_block_label(const int*icol, const char*name);

  ///
  /// tell splash about the label of a particle type
  ///
  /// \param[in] ityp  column of data
  /// \param[in] name  label for this column
  ///
  /// \note called by read_falcON_snapshot()
  void set_splash_particle_label(const int*ityp, const char*name);

  //
  //  1.2  routines provided here, but to be called from SPLASH
  //

  ///
  /// set debugging level
  ///
  /// \param[in] debug  debugging level
  ///
  /// \note  currently only debug=0,1,2,3 are distinguished.
  ///        0 means no output in case of an error,
  ///        1 means diagnostic output in case of error,
  ///        2 also prints errors from the HDF5 library
  ///        3 may print extra information (used in debugging this file)
  void set_falcON_debugging_level(const int*debug);

  ///
  /// opens a falcON HDF5 snapshot file
  ///
  /// \param[in]  filename       name of data file
  /// \param[out] ierr           non-zero if file couldn't be opened
  ///
  /// \note  precondition: none
  /// \note  postcondition: open_falcON_snapshot() can be called
  /// \note  We close any previously opened file (only one snapshot file can be
  ///        open at any time with this implementation), but warn if the
  ///        filename matches with the currently open file, if any.
  /// \note  Use num_falcON_snapshots() for the number of snapshots in the file
  void open_falcON_file(const char*filename, int*ierr);

  ///
  /// queries whether a file is open
  ///
  int falcON_file_is_open();

  ///
  /// closes the currently open (if any) falcON HDF5 snapshot file.
  ///
  /// \note  precondition: none
  /// \note  postcondition: falcON_file_is_open() will return false
  void close_falcON_file();

  ///
  /// queries if there is another snapshot present the currently open file
  ///
  /// \note  precondition: falcON_file_is_open()
  int have_falcON_snapshot(int*ierr);

  ///
  /// query the number of snapshots in the currently open file
  ///
  /// \note  precondition: falcON_file_is_open()
  /// \note  When using have_falcON_snapshot(), it may not be necessary to call
  ///        this function, see main() in PART 3 below for an example.
  /// \note  This function is trivial for falcON HDF5-based snapshot files,
  ///        but non-trivial for NEMO snapshot files, which we may support in
  ///        the future without changing this interface.
  int num_falcON_snapshots(int*ierr);

  ///
  /// opens the next falcON snapshot in the currently open file and
  /// reads its header
  ///
  /// \param[out] ntype  number of particle types
  /// \param[out] npart  number of particles per type
  /// \param[out] ncol   number of splash columns
  /// \param[out] dimX   number of spatial dimensions
  /// \param[out] dimV   number of velocity dimensions
  /// \param[out] time   simulation time of snapshot
  /// \param[out] hper   if hper[d]!=0, dimension d is periodic with |x|<hper
  /// 
  /// \note  precondition: have_falcON_snapshot()
  /// \note  postcondition: read_falcON_snapshot() can be called
  /// \note  particles types are (currently): gas,sink,std
  /// \note  closes any previously opened snasphots: this implemention only
  ///        allows one open snapshot at any time.
  /// \note  corresponds to 'read_XXX_header()' in systems where each snapshot
  ///        lives in its own file.
  void open_falcON_snapshot(int*ntype, int npart[max_num_types],
			    int*ncol, int*dimX, int*dimV, double*time,
			    double hper[3], int*ierr);

  ///
  /// read all data present in the snapshot
  ///
  /// \note  precondition: a snapshot has been opened via open_falcON_snapshot()
  /// \note  calls all functions declared in interface 1.1
  /// \note  generates splash column indices 0 to ncol-1
  /// \note  see also IMPORTANT NOTE below
  /// \note  tested in PART 3
  void read_falcON_snapshot(int*ierr);

  //
  //  IMPORTANT NOTE on falcON fields versus splash columns.
  //
  //  With the exception of field 'krnH', each scalar falcON field generates
  //  exactly one splash column. The column name is identical to the field name,
  //  except for field 'key', when the column name is 'id'.
  //  For gas particles, field 'krnH' becomes column 'h' with values reduced by
  //  a factor 2 to account for the convention in splash. For sink particles,
  //  field 'krnH' becomes column 'snkR', denoting the sink radius.
  //
  //  Vector fields, create 3 columns with names consisting of the first letter
  //  plus 'x', 'y', 'z', i.e. 'az' for the z-component of the acceleration.
  //  The only exceptions are the positions with column names 'x', 'y', 'z'.
  //
  //  Tensor fields, create 6 columns with names consisting of the capitalised
  //  first letter followed by 'xx', 'xy', 'xz', 'yy', 'yz', 'zz', i.e. 'Txy'
  //  for a component of the tidal field.
  //

} // extern "C"
////////////////////////////////////////////////////////////////////////////////
//
// PART 2  C++ implemention of interface 1.2
//
////////////////////////////////////////////////////////////////////////////////
#if !defined(__clang__) && defined(__GNUC__) && __GNUC__ >= 5
//
//  The HDF5 library may have been compiled with another (older) GNU C++ ABI
//  than the one used in the rest of the falcON.2 code. Here, we must use the
//  same C++ ABI as used with the HDF5 C++ library. This is controlled by the
//  following macro, see also
//
//     https://gcc.gnu.org/onlinedocs/libstdc++/manual/using_dual_abi.html
//
#  define _GLIBCXX_USE_CXX11_ABI 0
#endif
#if __cplusplus < 201103L
#  error requiring C++11
#endif
#include <iostream>
#include <cassert>
#include <type_traits>
#include <memory>
#include <vector>
#include <array>
#include <map>
#include <H5Cpp.h>
///
/// 2.1  auxiliary functionality
///      (symbols from the anonymous namespace have internal linkage
///
namespace {
  /// simple struct with data for each set of particles
  struct Particles {
    std::unique_ptr<H5::Group> group;        ///< HDF5 group
    std::string name;                        ///< name: 'sink', 'gas', 'std'
    std::size_t number;                      ///< number of particles
    // (I couldn't figure out to avoid this via aggregate initialisation)
    Particles(H5::Group*g, std::string const&t, std::size_t n)
      : group(g), name(t), number(n) {}
  };
  //
  using name_count_map = std::map<std::string,hsize_t>;
  using name_and_count = name_count_map::value_type;
  // static data
  int Debug = 0;                             ///< debug level (0 or not 0)
  int IndexSnap = 0;                         ///< index of next snapshot
  int NumSnap = 0;                           ///< number of snapshots in file
  int SplashCol = 0;                         ///< column for splash data
  int NumCol = 0;                            ///< number of splash columns
  name_count_map Fields;                     ///< field:dims map of all fields
  std::vector<Particles> Part;               ///< data for particle sets
  std::vector<double> Buffer;                ///< data buffer for input
  std::unique_ptr<H5::H5File> File;          ///< file
  std::string FileName;                      ///< name of currently open file
  const hsize_t k[1] = {3};
  const H5::ArrayType VectorType             ///< HDF5 data type for double[3]
  {H5::PredType::NATIVE_DOUBLE,1,k};
  // closes a splash column
  inline void close_column(const char*name)
  {
    if(SplashCol >= NumCol)
      throw "number of columns exceeds expected " + std::to_string(NumCol);
    set_splash_block_label(&SplashCol,name);
    SplashCol++;
  }
  // closes a splash column
  inline void close_column(std::string const&name)
  {
    close_column(name.c_str());
  }
  // closes currently open snapshot
  inline void close_snapshot()
  {
    SplashCol = 0;
    NumCol = 0;
    Part.clear();
    Fields.clear();
  }
  // select field for this particle type?
  // NOTE  Without selecting, we would potentially read meaningless data. This
  //       will not be necessary in future versions of falcON.
  inline bool select(std::string const&field, std::string const&ptype)
  {
    return
      (field=="spin" || field=="eabs" || field=="maxA")?
      (ptype=="sink") :
      (field=="snum" || field=="uin"  || field=="entr" || field=="dlKt" ||
       field=="dlKe" || field=="srho" || field=="alfa" || field=="divv" ||
       field=="dlht" || field=="vsig" || field=="fact" || field=="csnd" ||
       field=="pres" || field=="vort" || field=="dtdv" || field=="qmin" ||
       field=="delE" || field=="coll")?
      (ptype=="gas") :
      (field=="krnH" || field=="maxR") ?
      (ptype=="sink" || ptype=="gas") : true;
  }
  // try to read a component of a falcON field into a splash column
  void read_column(std::string const&field, const hsize_t comp,
		   const hsize_t dims)
  {
    bool read = false;
    // loop particle types
    for(int type=0; type!=int(Part.size()); ++type) {
      const auto&part = Part[std::size_t(type)];
      if(!select(field,part.name))
	continue;
      H5::DataSet data;
      try {
	data = part.group->openDataSet(field);
      } catch(...) {
	continue;  // field not present: continue with next particle type
      }
      // obtain size of data set
      auto space = data.getSpace();
      hsize_t count[2];
      auto rank = space.getSimpleExtentDims(count);
      assert(part.number==count[0]);
      const int ndat = int(part.number);
      // read component
      Buffer.resize(part.number);
      switch(dims) {
      case 1:
	// scalar field
	assert(rank==1);
	data.read(Buffer.data(),H5::PredType::NATIVE_DOUBLE,space,space);
	// reduce smoothing length by factor 2
	if(field=="krnH" && part.name=="gas")
	  for(auto&x:Buffer) x*=0.5;
	if(Debug>2) std::clog<<"reading "<<ndat<<" '"<<field<<"' for '"
			     <<part.name<<"' (type="<<type<<") into col="
			     <<SplashCol<<std::endl;
	read_falcON_data_into_splash(&SplashCol,&ndat,Buffer.data(),&type);
	read = true;
	// account for special case 'krnH'
	if(field=="krnH") {
	  if(part.name=="gas") {
	    close_column("h");
	    read = false;
	  } else if(part.name=="sink") {
	    close_column("snkR");
	    read = false;
	  }
	}
	break;
      default: {
	// component of vector/tensor field
	assert(rank==2);
	assert(dims==count[1]);
	count[1] = 1;
	const hsize_t start[2]={0,comp};
	space.selectHyperslab(H5S_SELECT_SET,count,start);
	data.read(Buffer.data(),H5::PredType::NATIVE_DOUBLE,
		  H5::DataSpace{1,count},space);
	if(Debug>2) std::clog<<"reading "<<ndat<<" '"<<field<<'['<<comp
			     <<"]' for '"<<part.name<<"' (type="<<type
			     <<") into col="<<SplashCol<<std::endl;
	read_falcON_data_into_splash(&SplashCol,&ndat,Buffer.data(),&type);
	read = true;
      }
      }
    }
    // if any data have been read: set column label
    if(read) {
      std::string name;
      if(field=="pos")
	name = comp==0? "x" : comp==1? "y" : "z";
      else if(dims==3) {
	name = field[0];
	name+= comp==0? "x" : comp==1? "y" : "z";
      } else if(dims==6) {
	name = static_cast<char>(std::toupper(field[0]));
	name+= comp==0? "xx": comp==1? "xy": comp==2? "xz":
	       comp==3? "yy": comp==4? "yz": "zz";
      } else
	name=field;
      close_column(name);
    }
  }
  // read a field
  inline void read_field(name_and_count const&field)
  {
    for(hsize_t comp=0; comp!=field.second; ++comp)
      read_column(field.first, comp, field.second);
  }
  // try to read a falcON field
  inline void read_field(std::string const&field)
  {
    const auto field_iter = Fields.find(field);
    if(field_iter==Fields.end())
      std::clog<<"WARNING: falcON field '"<<field<<"' not present in file"
	       <<std::endl;
    else
      read_field(*field_iter);
  }
}
//
//  2.2  implement routines declared in 1.2
//
void set_falcON_debugging_level(const int*d)
{
  Debug = *d;
}
//
void open_falcON_file(const char*file, int*ierr)
{
  *ierr = 1;
  IndexSnap = 0;
  close_snapshot();
  // suppress error messages from the HDF5 C library
  if(Debug<2) H5::Exception::dontPrint();
  if(file==nullptr) {
    if(Debug) std::clog<<"open_falcON_file(): empty filename"<<std::endl;
    return close_falcON_file();
  }
  if(FileName==file)
    std::clog<<"open_falcON_file('"<<file<<"'): WARNING: "
	     <<"a falcON hdf5 snapshot file of the same name is already open, "
	     <<"re-opening it will effectively rewind to the first snapshot."
	     <<std::endl;
  // open HDF5 file
  FileName = file;
  try {
    File.reset(new H5::H5File(file,H5F_ACC_RDONLY));
  } catch(H5::Exception const&exc) {
    if(Debug) std::clog<<"open_falcON_file(): "
		       <<"cannot open HDF5 file '"<<file<<"' (HDF5 says \""
		       <<exc.getDetailMsg()<<"\")"<<std::endl;
    return close_falcON_file();
  }
  // open attribute 'falcON' to establish that this is a falcON snapshot file
  try {
    auto attr = File->openAttribute("falcON");
  } catch(...) {
    if(Debug) std::clog<<"open_falcON_file(): file '"<<file
		       <<"' is an HDF5 file, but not a falcON snapshot file"
		       <<std::endl;
    return close_falcON_file();
  }
  // read number of snapshots
  try {
    auto attr = File->openAttribute("num_snapshots");
    attr.read(H5::PredType::NATIVE_UINT32,&NumSnap);
  } catch(H5::Exception const&exc) {
    // should never happen
    if(Debug) std::clog<<"open_falcON_file('"<<file<<"'): failed to "
		       <<"read num_snapshots (HDF5 says \""
		       <<exc.getDetailMsg()<<"\")"<<std::endl;
    return close_falcON_file();
  }
  *ierr = 0;
}
//
int falcON_file_is_open()
{
  return bool(File);
}
//
void close_falcON_file()
{
  close_snapshot();
  IndexSnap = NumSnap = 0;
  File.reset();
  FileName.clear();
}
//
int have_falcON_snapshot(int*ierr)
{
  *ierr = !File;
  if(*ierr && Debug)
    std::clog<<"have_falcON_snapshot(): no snapshot file open"<<std::endl;
  return IndexSnap < NumSnap;
}
//
int num_falcON_snapshots(int*ierr)
{
  *ierr = !File;
  if(*ierr && Debug)
    std::clog<<"num_falcON_snapshots(): no snapshot file open"<<std::endl;
  return NumSnap;
}
//
void open_falcON_snapshot(int*ntyp, int npart[max_num_types],
			  int*ncol, int*dimX, int*dimV,
			  double*time, double hper[3], int*ierr)
{
  *ierr = 1;
  *ntyp = 0;
  *ncol = 0;
  *dimX = 3;
  *dimV = 3;
  close_snapshot();
  // check preconditions
  if(!File) {
    if(Debug) std::clog<<"open_falcON_snapshot(): "
		       <<"no snapshot file open"<<std::endl;
    return;
  }
  if(IndexSnap >= NumSnap) {
    if(Debug) std::clog<<"open_falcON_snapshot(): no more than "
		       <<NumSnap<<" snapshots in file"<<std::endl;
    return;
  }
  try {
    // open snapshot group
    auto name = "snapshot" + std::to_string(IndexSnap++);
    auto snap = File->openGroup(name);
    // read time and hper
    auto attr = snap.openAttribute("time");
    attr.read(H5::PredType::NATIVE_DOUBLE,time);
    attr = snap.openAttribute("hper");
    attr.read(VectorType,hper);
    // read npart[] and open particle sets
    std::array<std::string,3> types = {{"sink","gas","std"}};
    for(const auto&type:types) {
      name = "N";
      name+= type;
      unsigned number;
      try {
	attr = snap.openAttribute(name);
	attr.read(H5::PredType::NATIVE_UINT32,&number);
      } catch(...) {
	// should never go here
	number = 0;
      }
      if(number) {
	// open HDF5 group for particles
	const auto part = snap.openGroup(type);
	Part.emplace_back(new H5::Group(part),type,number);
	npart[(*ntyp)++] = int(number);
	// collect (field:dims) pairs and count columns
	for(hsize_t fld=0; fld!=part.getNumObjs(); ++fld) {
	  const auto field = part.getObjnameByIdx(fld);
	  const auto have_field = Fields.count(field);
	  if(have_field) {
	    // field 'krnH' generates an extra column for each particle type
	    if(field=="krnH") (*ncol)++;
	  } else {
	    hsize_t count[2];
	    const auto dims = 1==part.openDataSet(field).getSpace().
	      getSimpleExtentDims(count)? 1:count[1];
	    Fields[field] = dims;
	    (*ncol)+= int(dims);
	  }
	  if(Debug>2) std::clog<<"counting columns: type="<<type
			       <<" fld="<<fld<<"='"<<field<<"' cols->"
			       <<(*ncol)<<std::endl;
	}
      }
    }
  } catch(H5::Exception const&exc) {
    // should never go here
    if(Debug) std::clog<<"open_falcON_snapshot(): HDF5 error: \""
		       <<exc.getDetailMsg()<<'\"'<<std::endl;
    *ntyp = 0;
    *ncol = 0;
    return close_snapshot();
  }
  NumCol = *ncol;
  *ierr = 0;
}
//
void read_falcON_snapshot(int*ierr)
{
  *ierr = 0;
  std::vector<std::string> PreferredOrder = 
    {"pos", "key", "vel", "acc", "mass", "pot", "pex", "rung",
     "krnH", "srho", "uin", "entr", "divv", "dlKt", "dlKe", "alfa"};
  try {
    auto FieldsToRead=Fields;
    // read some fields in preferred order
    for(const auto&field:PreferredOrder)
      if(FieldsToRead.erase(field))
	read_field(field);
    // read remaining fields
    for(const auto&field:FieldsToRead)
      read_field(field);
  } catch(H5::Exception const&exc) {
    // catch any HDF5 error
    if(Debug) std::clog<<"read_falcON_snapshot(): HDF5 error: \""
		       <<exc.getDetailMsg()<<'\"'<<std::endl;
    return;
  } catch(std::string const&exc) {
    // catch error from column overflow
    if(Debug)
      std::clog<<exc<<" -- this either means that read_falcON_snapshot() "
	       <<"has been called more than once, or that "
	       <<"there is an error in file read_falcON_utils.cc"<<std::endl;
    return;
  }
  // set particle type names
  for(int type=0; type!=int(Part.size()); ++type)
    set_splash_particle_label(&type, Part[std::size_t(type)].name.c_str());
  *ierr = 0;
}
#if(0) // untested
//
void query_falcON_fields(int*nfld, char*flds, int*ierr)
{
  *ierr = 1;
  *nfld = 0;
  // check precondition
  if(Part.empty()) {
    if(Debug) std::clog<<"query_falcON_fields(): no snapshot open"<<std::endl;
    return;
  }
  // set particle type names
  for(int type=0; type!=int(Part.size()); ++type)
    set_splash_particle_label(&type, Part[std::size_t(type)].name.c_str());
  // copy fields
  for(const auto&field:Fields) {
    assert(field.first.size() < 5);
    assert(*nfld < max_num_cols);
    strcpy(flds,field.first.c_str());
    flds+= 5;
    (*nfld)++;
  }
  *ierr = 0;
}
//
void read_falcON_field(const char*field, int*ierr)
{
  *ierr = 1;
  if(field==nullptr) {
    if(Debug)
      std::clog<<"read_falcON_field(): empty field"<<std::endl;
    return;
  }
  try {
    read_field(field);
  } catch(H5::Exception const&exc) {
    // catch any HDF5 error
    if(Debug) std::clog<<"read_falcON_field(): HDF5 error: \""
		       <<exc.getDetailMsg()<<'\"'<<std::endl;
    return;
  } catch(std::string const&exc) {
    // catch error from column overflow
    if(Debug)
      std::clog<<exc<<std::endl;
    return;
  }
  *ierr = 0;
}
#endif
#ifdef TestMain
////////////////////////////////////////////////////////////////////////////////
//
// PART 3  test above implementation
//
////////////////////////////////////////////////////////////////////////////////
#include <iomanip>
#include <cstdlib>
#include <algorithm>
//
// 3.1  implement interface 1.1 using C++
//
namespace {
  struct data_per_particle_type {
    std::string name;
    std::size_t number;
    std::vector< std::vector<double> > columns;
  };
  std::vector<data_per_particle_type> particle_data;
  std::vector<std::string> column_labels;
}
//
void set_splash_block_label(const int*icol, const char*name)
{
  column_labels.at(std::size_t(*icol)) = name;
}
//
void set_splash_particle_label(const int*type, const char*name)
{
  particle_data.at(std::size_t(*type)).name = name;
}
//
void read_falcON_data_into_splash(const int*icol, const int*ndat,
				  const double*from, const int*type)
{
  auto&part = particle_data.at(std::size_t(*type));
  assert(*ndat > 0);
  assert(part.number == std::size_t(*ndat));
  part.columns.at(std::size_t(*icol)).resize(std::size_t(*ndat));
  std::copy(from,from+*ndat,part.columns[std::size_t(*icol)].begin());
}
//
// 3.2  an executable that dumps a falcON file
//
int main(int argc, const char**argv)
{
  // obtain program parameters
  if(argc<2) {
    std::clog<<"usage: '"<<argv[0]
             <<" input_file [lines_per_type [debug_level]]'"<<std::endl;
    return 0;
  }
  const auto input = argv[1];
  const auto lines = argc<3? 10 : std::atoi(argv[2]);
  const auto debug = argc<4?  0 : std::atoi(argv[3]);
  if(lines<=0) {
    std::clog<<"ERROR: lines="<<lines<<"<=0"<<std::endl;
    return 1;
  }
  if(debug<0) {
    std::clog<<"ERROR: debug="<<debug<<"<0"<<std::endl;
    return 1;
  }
  // set debug level
  set_falcON_debugging_level(&debug);
  // open falcON snapshot file
  int ierr;
  open_falcON_file(input,&ierr);
  if(ierr) {
    std::clog<<"ERROR: open_falcON_file() failed"<<std::endl;
    return 1;
  }
  if(!have_falcON_snapshot(&ierr)) {
    std::clog<<"WARNING: have_falcON_snapshot()==false after "
	     <<"open_falcON_file()"<<std::endl;
    return 1;
  }
  // loop snapshots
  while(have_falcON_snapshot(&ierr)) {
    // open snapshot and read header
    int dimX, dimV, ntype, ncol, npart[max_num_types];
    double time, hper[3];
    open_falcON_snapshot(&ntype,npart,&ncol,&dimX,&dimV,&time,hper,&ierr);
    if(ierr) {
      std::clog<<"ERROR: open_falcON_snapshot() failed"<<std::endl;
      return 1;
    }
    if(ntype<=0) {
      std::clog<<"WARNING: open_falcON_snapshot() returned ntype<=0"
	       <<std::endl;
      return 1;
    }
    if(ncol<=0) {
      std::clog<<"WARNING: open_falcON_snapshot() returned ncol<=0"
	       <<std::endl;
      return 1;
    }
    if(debug)
      std::clog<<"DEBUG: open_falcON_snapshot(): ntype="<<ntype
	       <<" ncol="<<ncol<<std::endl;
    // prepare particle data for reading
    particle_data.resize(std::size_t(ntype));
    column_labels.resize(std::size_t(ncol));
    for(int type=0; type!=ntype; ++type) {
      auto&part = particle_data.at(std::size_t(type));
      assert(npart[type]>0);
      part.number = std::size_t(npart[type]);
      part.columns.resize(std::size_t(ncol));
    }
    // read particle data
    if(debug)
      std::clog<<"DEBUG: now calling read_falcON_snapshot()"
	       <<std::endl;
    read_falcON_snapshot(&ierr);
    if(ierr) {
      std::clog<<"ERROR: read_falcON_snapshot() failed"<<std::endl;
      return 1;
    }
    // print header for columns
    for(int c=0; c!=ncol; ++c)
      std::cout<<"=============";
    std::cout<<'\n';
    std::cout<<" time="<<time
             <<" hper=["<<hper[0]<<' '<<hper[1]<<' '<<hper[2]<<']';
    for(int type=0; type!=ntype; ++type)
      std::cout<<" n_"<<Part[std::size_t(type)].name<<'='<<npart[type];
    std::cout<<'\n';
    for(int c=0; c!=ncol; ++c)
      std::cout<<"-------------";
    std::cout<<'\n';
    for(int c=0; c!=ncol; ++c)
      std::cout<<"       col "<<std::setw(2)<<c;
    std::cout<<'\n';
    for(int c=0; c!=ncol; ++c)
      std::cout<<' '<<std::setw(12)<<column_labels.at(std::size_t(c));
    std::cout<<'\n';
    // loop particle types and print columns
    for(int type=0; type!=ntype; ++type) {
      const auto&part = particle_data.at(std::size_t(type));
      const auto n = lines<npart[type]? lines : npart[type];
      for(int c=0; c!=ncol; ++c)
	std::cout<<"-------------";
      std::cout<<'\n';
      std::cout<<" '"<<part.name<<"' particles:\n";
      for(int c=0; c!=ncol; ++c)
	std::cout<<"-------------";
      std::cout<<'\n';
      for(int i=0; i<n; ++i) {
	for(int c=0; c!=ncol; ++c) {
	  const auto&column = part.columns.at(std::size_t(c));
	  if(column.empty()) std::cout<<"            -";
	  else std::cout<<' '<<std::setw(12)<<column.at(std::size_t(i));
	}
	std::cout<<'\n';
      }
      if(n<npart[type]) {
	for(int c=0; c!=ncol; ++c)
          std::cout<<"          ...";
	std::cout<<'\n';
      }
      std::cout<<std::flush;
    }
  }
  // close snapshot file
  close_falcON_file();
}
#endif
