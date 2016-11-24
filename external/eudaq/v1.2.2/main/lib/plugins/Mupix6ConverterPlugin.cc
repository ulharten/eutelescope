/**
 * an eudaq converter plugin for the MuPix6 readout
 *
 * @author      Heiko Augustin <augustin@physi.uni-heidelberg.de>
 * @date        2015-02
 *
 */

#ifndef WIN32
#include <inttypes.h> /* uint32_t */
#endif
#include <memory>
#include <vector>

#include "eudaq/DataConverterPlugin.hh"
#include "eudaq/Logger.hh"
#include "eudaq/Mupix.hh"
#include "eudaq/RawDataEvent.hh"
#include "eudaq/StandardEvent.hh"
#include "eudaq/Utils.hh"
#include "eudaq/telescope_frame.hh"

#if USE_LCIO
#include "IMPL/LCEventImpl.h"
#include "IMPL/TrackerRawDataImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "lcio.h"
#endif

#if USE_EUTELESCOPE
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTakiDetector.h"
#include "EUTelSetupDescription.h"
#include "EUTelEventImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelTrackerDataTriggerInterfacer.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelMuPixel.h"
#include "EUTelExternalTrigger.h"
#endif

using namespace std;

static const char *     MUPIX_EVENT_TYPE = "MUPIX7";
static const unsigned   MUPIX_SENSOR_ID = 71;
static const char *     MUPIX_SENSOR_TYPE = "MUPIX7";
static const unsigned   MUPIX_TYPE = 7;
static const unsigned   MUPIX_SENSOR_NUM_COLS = 40;
static const unsigned   MUPIX_SENSOR_NUM_ROWS = 32;
static const unsigned   MUPIX_SENSOR_BINARY_SIGNAL = 1;
static const char *     MUPIX_COLLECTION_NAME = "zsdata_mupix7";
static const char *     TRIGGER_COLLECTION_NAME = "eudet_triggers";
static const char *     TOT_COLLECTION_NAME = "eudet_tots";


template <typename T> int sgn(T val) {
  return  (T(0) < val)  - ( val < T(0) ); 
}


namespace eudaq {

class Mupix6ConverterPlugin : public DataConverterPlugin {
public:

    virtual void Initialize(const eudaq::Event &, eudaq::Configuration const &) {};
    virtual unsigned GetTriggerID(const eudaq::Event &) const;
    virtual bool GetStandardSubEvent(StandardEvent &, const eudaq::Event &) const;

#if USE_LCIO && USE_EUTELESCOPE
    // run header conversion is not implemented since it seems to be not
     // really used anymore by the EUTelescope framework
  virtual bool GetLCIOSubEvent(lcio::LCEvent &, const eudaq::Event *,  const eudaq::Event *) const;
  virtual bool GetLCIOSubEvent(lcio::LCEvent &, const eudaq::Event *,  const eudaq::Event *, const eudaq::Event *) const;
#endif

private:
    // The constructor can be private since only one static instance is created
    // The DataConverterPlugin constructor must be passed the event type
    // in order to register this converter for the corresponding conversions
    // Member variables should also be initialized to default values here.
    Mupix6ConverterPlugin() : DataConverterPlugin(MUPIX_EVENT_TYPE) {}

    // static instance as mentioned above
    static Mupix6ConverterPlugin m_instance;
};

/** create the singleton. registers the converter with the plugin-manager */
Mupix6ConverterPlugin Mupix6ConverterPlugin::m_instance;

/** return the trigger id (as provided by the TLU).
 *
 * returns (unsigned)(-1) if it can not be retrieved.
 */
unsigned Mupix6ConverterPlugin::GetTriggerID(const Event & ev) const
{
    if (const RawDataEvent * rev = dynamic_cast<const RawDataEvent *> (&ev)){
    // trigger id is not defined for special events
    if (rev->IsBORE() || rev->IsEORE()) return (unsigned)(-1);
    if (rev->GetID(rev->NumBlocks()-1)==unsigned(-1)) return (unsigned)(-1);
    else return rev->GetID(rev->NumBlocks()-1);
    }
}

/** convert the data in the RawDataEvent into a StandardEvent
 *
 * should return false on failure but this does not seem to have any effect.
 */
bool Mupix6ConverterPlugin::GetStandardSubEvent(
    StandardEvent & dest,
    const Event & ev) const
{
    
    uint col=0, row=0;
    uint64_t nhits=0;
    
    if (const RawDataEvent * source = dynamic_cast<const RawDataEvent *> (&ev)){
    TelescopeFrame data=TelescopeFrame();
        
    // beginning-of-run / end-of-run should not be converted
    if (source->IsBORE()) {
        // this should never happen. BORE event should be handled
        // specially by the initialize function
        EUDAQ_ERROR("got BORE during conversion");
        return true;
    } else if (source->IsEORE()) {
        // i'm not sure if this can happen
        EUDAQ_WARN("got EORE during conversion");
        return true;
    }

    StandardPlane plane(MUPIX_SENSOR_ID, MUPIX_EVENT_TYPE, MUPIX_SENSOR_TYPE);
    
    for( uint i=0 ; i<source->NumBlocks() ; i++){
      
        data.from_bytes(source->GetBlock(i));
        nhits+=data.num_hits();
        data.clear();
    }
    
    plane.SetSizeZS(MUPIX_SENSOR_NUM_COLS, MUPIX_SENSOR_NUM_ROWS,nhits);
    plane.SetTLUEvent(source->GetID(source->NumBlocks()-1));
    //std::cout<<"MupixConverter "<<source->GetID(source->NumBlocks()-1)<<std::endl;
    nhits=0;
    for( uint i=0 ; i<source->NumBlocks() ; i++){
      
	data.from_bytes(source->GetBlock(i));
	
	for (unsigned j = 0; j < data.num_hits() ; ++j){
            
            col=data.get_hit(j,MUPIX_TYPE).column();
            row=data.get_hit(j,MUPIX_TYPE).row();
	    if(!(col==0&&row==0)||col>31 ||row>39) plane.SetPixel(nhits+j,row ,col ,MUPIX_SENSOR_BINARY_SIGNAL,0,0);
	
	}
	nhits+=data.num_hits();
	data.clear();
    }

    // Furthermore we put the timestamp into the signal amplitude
    if(source->GetID(source->NumBlocks()-1)>100) dest.AddPlane(plane);
    else{
        StandardPlane plane2(MUPIX_SENSOR_ID, MUPIX_EVENT_TYPE, MUPIX_SENSOR_TYPE);
        plane2.SetSizeZS(MUPIX_SENSOR_NUM_COLS, MUPIX_SENSOR_NUM_ROWS,0);
        plane2.SetTLUEvent(source->GetID(source->NumBlocks()-1));
        dest.AddPlane(plane2);
    }

    return true;
    }
}

#if USE_LCIO && USE_EUTELESCOPE 

bool Mupix6ConverterPlugin::GetLCIOSubEvent(
    lcio::LCEvent & dest,
    const eudaq::Event * ev,
    const eudaq::Event * next) const
{
    using lcio::CellIDEncoder;
    using lcio::LCEvent;
    using lcio::LCCollectionVec;
    using lcio::TrackerDataImpl;
    using eutelescope::EUTELESCOPE;
    using eutelescope::EUTelTrackerDataInterfacerImpl;
    using eutelescope::EUTelMuPixel;
    using eutelescope::EUTelExternalTrigger;
    using eutelescope::EUTelTrackerDataTriggerInterfacer;
        
    
    // beginning-of-run / end-of-run should not be converted
    if (ev->IsBORE()) {
        // this should never happen. BORE event should be handled
        // specially by the initialize function
        EUDAQ_ERROR("got BORE during lcio conversion");
        return true;
    } else if (ev->IsEORE()) {
        // i'm not sure if this can happen
        EUDAQ_WARN("got EORE during lcio conversion");
        return true;
    }

    // lcio takes ownership over the different data objects when they are
    // added to the event. secure them w/ auto_ptr so they get deleted
    // automatically when something breaks.
    // NOTE to future self: use unique_ptr since auto_ptr is deprecated
    bool mupix_collection_exists = false;
    bool trigger_collection_exists = false;
    bool tot_collection_exists = false;
    LCCollectionVec * mupix_collection;
    LCCollectionVec * trigger_collection;
    LCCollectionVec * tot_collection;
    std::auto_ptr<TrackerDataImpl> mupix_frame;
    std::auto_ptr<TrackerDataImpl> trigger_frame;
    std::auto_ptr<TrackerDataImpl> tot_frame;
    std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelMuPixel> > pixels;
    std::auto_ptr<EUTelTrackerDataTriggerInterfacer> triggers;
    std::auto_ptr<EUTelTrackerDataTriggerInterfacer> tots;

    // get the lcio output collections
    try {
        mupix_collection = static_cast<LCCollectionVec *>(
            dest.getCollection(MUPIX_COLLECTION_NAME));
        mupix_collection_exists = true;
    } catch(lcio::DataNotAvailableException & e) {
        mupix_collection = new LCCollectionVec(lcio::LCIO::TRACKERDATA);
        mupix_collection_exists = false;
    }
    try {
        trigger_collection = static_cast<LCCollectionVec *>(
            dest.getCollection(TRIGGER_COLLECTION_NAME));
        trigger_collection_exists = true;
    } catch(lcio::DataNotAvailableException & e) {
        trigger_collection = new LCCollectionVec(lcio::LCIO::TRACKERDATA);
        trigger_collection_exists = false;
    }
    try {
        tot_collection = static_cast<LCCollectionVec *>(
            dest.getCollection(TOT_COLLECTION_NAME));
        tot_collection_exists = true;
    } catch(lcio::DataNotAvailableException & e) {
        tot_collection = new LCCollectionVec(lcio::LCIO::TRACKERDATA);
        tot_collection_exists = false;
    }
    

    // a cell id identifies from which detector component a specific data
    // collection originates.
    // here it is used to encode the telescope / dut plane and the type of
    // data, e.g. zero-suppressed or raw, that is stored in the collection
    CellIDEncoder<TrackerDataImpl>
        cell_encoder(EUTELESCOPE::ZSDATADEFAULTENCODING, mupix_collection);
    if (MUPIX_SENSOR_ID == 601) {
      cell_encoder["sensorID"] = 61;
    }
    else if (MUPIX_SENSOR_ID == 701) {
      cell_encoder["sensorID"] = 71;
    }
    else {
      cell_encoder["sensorID"] = MUPIX_SENSOR_ID;
    }
    cell_encoder["sparsePixelType"] = eutelescope::kEUTelMuPixel;
    
    CellIDEncoder<TrackerDataImpl>
      cell_encoder_trigger(EUTELESCOPE::ZSDATADEFAULTENCODING, trigger_collection);
    cell_encoder_trigger["sensorID"] = 1;

    CellIDEncoder<TrackerDataImpl>
      cell_encoder_tot(EUTELESCOPE::ZSDATADEFAULTENCODING, tot_collection);
    cell_encoder_tot["sensorID"] = 1;
    
    // the lcio objects that store the hit and trigger data for a single readout frame
    mupix_frame.reset(new TrackerDataImpl);
    cell_encoder.setCellID(mupix_frame.get());
    trigger_frame.reset(new TrackerDataImpl);
    cell_encoder_trigger.setCellID(trigger_frame.get());
    tot_frame.reset(new TrackerDataImpl);
    cell_encoder_trigger.setCellID(tot_frame.get());

    // a convenience object that encodes the sparse pixel data and external triggers into an
    // eutelescope-specific format and stores them in the given readout frames
    pixels.reset( new EUTelTrackerDataInterfacerImpl<EUTelMuPixel>( mupix_frame.get() ) );
    triggers.reset( new EUTelTrackerDataTriggerInterfacer ( trigger_frame.get() ) );
    tots.reset( new EUTelTrackerDataTriggerInterfacer ( tot_frame.get() ) );

    if (const RawDataEvent * source = dynamic_cast<const RawDataEvent *> (ev)){
      if (const RawDataEvent * next_source = dynamic_cast<const RawDataEvent *> (next)){
	
	TelescopeFrame data=TelescopeFrame();
	//cout << "current: " << source->NumBlocks() << " next: " << next_source->NumBlocks() << endl;
	for ( uint i = 0; i < source->NumBlocks(); ++i ) {
	  data.from_bytes( source->GetBlock( i ) ); 
	  
	  /* Extract TLU trigger and other triggers
	     timestamp: RawTrigger::timestamp() returns the absolute timestamp, not the difference
	     between frame and trigger,
	     tag: 0x1 for TLU trigger, 0xBA for "normal" triggers
	     
	     Create ExternalTrigger object which stores timestamps (64 bits) and labels (16 bits)
	     of external triggers,
	     add to triggers LCIO data object
	  */
	  RawTrigger trig;
//	  cout << "size of triggers: " << data.num_triggers() << endl;
	  for ( uint j = 0; j < data.num_triggers(); j++) {
	    trig = data.get_trigger(j);
	    EUTelExternalTrigger eutel_trigger(
					       trig.timestamp(),
					       trig.tag()
					       );
	    triggers->addExternalTrigger( &eutel_trigger );
	  }
	  
	  /* Extract ToT information
	     Create ExternalTrigger object which stores timestamp (48 bits) plus length of ToT (8 bits) 
	     and a  label to identify as ToT (0x2),
	     add to triggers LCIO data object
	  */
	  RawTimeOverThreshold tot;
	  for ( uint j = 0; j < data.num_tots(); j++) {
	    tot = data.get_tot(j);
	    EUTelExternalTrigger eutel_tot(
					   (uint64_t)( (tot.timestamp() & 0xFFFFFFFFFFFF)<<8 | (tot.length() & 0xFF ) ),
					   0x2
	  				   );
	    tots->addExternalTrigger( &eutel_tot );
	  }
	  
	  /* 
	     Create MuPixel object which also stores the hit and frame timestamp,
	     add to pixels LCIO data object
	  */
	  for ( uint j = 0; j < data.num_hits(); ++j ) {
	    // Caution: columns and rows swapped for EUTelescope analysis
	    // of March 2016 DESY data
	    // since 90 degree rotation is not supported by the framework
	    EUTelMuPixel pixel(
			       //data.get_hit(j,MUPIX_TYPE).column(),
			       //data.get_hit(j,MUPIX_TYPE).row(),
			       data.get_hit(j,MUPIX_TYPE).row(),
			       data.get_hit(j,MUPIX_TYPE).column(),
			       1, // binary
			       0,
			       data.get_hit(j,MUPIX_TYPE).timestamp_raw(),  // hit time stamp (8 bit)
			       data.timestamp() & 0xFFFFFFFF  // frame time stamp, 64 bit!! LCIO can only store floats
			       );
	    if ( data.get_hit(j,MUPIX_TYPE).row() < MUPIX_SENSOR_NUM_COLS && data.get_hit(j,MUPIX_TYPE).column() < MUPIX_SENSOR_NUM_ROWS ) {
	      pixels->addSparsePixel( &pixel );
	    }
	    else 
	      cout << "WARNING: col = " << data.get_hit(j,MUPIX_TYPE).row() << ", row = " << data.get_hit(j,MUPIX_TYPE).column() << endl;
	    
	  }
	  data.clear();
	}
 	
	for ( uint i = 0; i < next_source->NumBlocks(); ++i ) {
	  data.from_bytes( next_source->GetBlock( i ) );

	  /* Extract TLU trigger and other triggers
	     timestamp: RawTrigger::timestamp() returns the absolute timestamp, not the difference
	     between frame and trigger,
	     tag: 0x1 for TLU trigger, 0xBA for "normal" triggers

	     Create ExternalTrigger object which stores timestamps (64 bits) and labels (16 bits)
	     of external triggers,
	     add to triggers LCIO data object
	  */
	  RawTrigger trig;
	  for ( uint j = 0; j < data.num_triggers(); j++) {
	    trig = data.get_trigger(j);
	    EUTelExternalTrigger eutel_trigger(
					       trig.timestamp(),
					       trig.tag()
					       );
	    triggers->addExternalTrigger( &eutel_trigger );
	  }

	  /* Extract ToT information
	     Create ExternalTrigger object which stores timestamp (48 bits) plus length of ToT (8 bits)
	     and a  label to identify as ToT (0x2),
	     add to triggers LCIO data object
	  */
	  RawTimeOverThreshold tot;
	  for ( uint j = 0; j < data.num_tots(); j++) {
	    tot = data.get_tot(j);
	    EUTelExternalTrigger eutel_tot(
	  				   (uint64_t)( (tot.timestamp() & 0xFFFFFFFFFFFF)<<8 | (tot.length() & 0xFF ) ),
					   0x2
	  				   );
	    tots->addExternalTrigger( &eutel_tot );
	  }

	  /*
	     Create MuPixel object which also stores the hit and frame timestamp,
	     add to pixels LCIO data object
	  */
	  for ( uint j = 0; j < data.num_hits(); ++j ) {
	    // Caution: columns and rows swapped for EUTelescope analysis
	    // since 90 degree rotation is not supported by the framework
	    EUTelMuPixel pixel(
			       //data.get_hit(j,MUPIX_TYPE).column(),
			       //data.get_hit(j,MUPIX_TYPE).row(),
			       data.get_hit(j,MUPIX_TYPE).row(),
			       data.get_hit(j,MUPIX_TYPE).column(),
			       1, // binary
			       0,
			       data.get_hit(j,MUPIX_TYPE).timestamp_raw(),  // hit time stamp (8 bit)
			       data.timestamp() & 0xFFFFFFFF  // frame time stamp, 64 bit!! LCIO can only store floats
			       );
	    if ( data.get_hit(j,MUPIX_TYPE).row() < MUPIX_SENSOR_NUM_COLS && data.get_hit(j,MUPIX_TYPE).column() < MUPIX_SENSOR_NUM_ROWS ) {
	      pixels->addSparsePixel( &pixel );
	    }
	    else
	      cout << "WARNING: col = " << data.get_hit(j,MUPIX_TYPE).row() << ", row = " << data.get_hit(j,MUPIX_TYPE).column() << endl;

	  }
	  data.clear();
	}
	
	// hand over ownership over the readout frame to the lcio collection
	mupix_collection->push_back( mupix_frame.release() );
	trigger_collection->push_back( trigger_frame.release() );
	tot_collection->push_back( tot_frame.release() );
	
	if (!mupix_collection_exists && (mupix_collection->size() != 0)) {
	  dest.addCollection(mupix_collection, MUPIX_COLLECTION_NAME);
	}
	else {
	  cout << "FAILED to convert Mupix event to RawDataEvent" << endl;
	}
	if (!trigger_collection_exists && (trigger_collection->size() != 0)) {
	  dest.addCollection(trigger_collection, TRIGGER_COLLECTION_NAME);
	}
	else {
	  cout << "FAILED to convert trigger event to RawDataEvent" << endl;
	}
	if (!tot_collection_exists && (tot_collection->size() != 0)) {
	  dest.addCollection(tot_collection, TOT_COLLECTION_NAME);
	}
	else {
	  cout << "FAILED to convert tot event to RawDataEvent" << endl;
	}
      }
    }
    
    return true;
}



// use three consecutive mupix readout cycles (searching for higher efficiency)

bool Mupix6ConverterPlugin::GetLCIOSubEvent(
    lcio::LCEvent & dest,
    const eudaq::Event * ev,
    const eudaq::Event * next,
	const eudaq::Event * nextNext) const
{
    using lcio::CellIDEncoder;
    using lcio::LCEvent;
    using lcio::LCCollectionVec;
    using lcio::TrackerDataImpl;
    using eutelescope::EUTELESCOPE;
    using eutelescope::EUTelTrackerDataInterfacerImpl;
    using eutelescope::EUTelMuPixel;
    using eutelescope::EUTelExternalTrigger;
    using eutelescope::EUTelTrackerDataTriggerInterfacer;


    // beginning-of-run / end-of-run should not be converted
    if (ev->IsBORE()) {
        // this should never happen. BORE event should be handled
        // specially by the initialize function
        EUDAQ_ERROR("got BORE during lcio conversion");
        return true;
    } else if (ev->IsEORE()) {
        // i'm not sure if this can happen
        EUDAQ_WARN("got EORE during lcio conversion");
        return true;
    }

    // lcio takes ownership over the different data objects when they are
    // added to the event. secure them w/ auto_ptr so they get deleted
    // automatically when something breaks.
    // NOTE to future self: use unique_ptr since auto_ptr is deprecated
    bool mupix_collection_exists = false;
    bool trigger_collection_exists = false;
    bool tot_collection_exists = false;
    LCCollectionVec * mupix_collection;
    LCCollectionVec * trigger_collection;
    LCCollectionVec * tot_collection;
    std::auto_ptr<TrackerDataImpl> mupix_frame;
    std::auto_ptr<TrackerDataImpl> trigger_frame;
    std::auto_ptr<TrackerDataImpl> tot_frame;
    std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelMuPixel> > pixels;
    std::auto_ptr<EUTelTrackerDataTriggerInterfacer> triggers;
    std::auto_ptr<EUTelTrackerDataTriggerInterfacer> tots;

    // get the lcio output collections
    try {
        mupix_collection = static_cast<LCCollectionVec *>(
            dest.getCollection(MUPIX_COLLECTION_NAME));
        mupix_collection_exists = true;
    } catch(lcio::DataNotAvailableException & e) {
        mupix_collection = new LCCollectionVec(lcio::LCIO::TRACKERDATA);
        mupix_collection_exists = false;
    }
    try {
        trigger_collection = static_cast<LCCollectionVec *>(
            dest.getCollection(TRIGGER_COLLECTION_NAME));
        trigger_collection_exists = true;
    } catch(lcio::DataNotAvailableException & e) {
        trigger_collection = new LCCollectionVec(lcio::LCIO::TRACKERDATA);
        trigger_collection_exists = false;
    }
    try {
        tot_collection = static_cast<LCCollectionVec *>(
            dest.getCollection(TOT_COLLECTION_NAME));
        tot_collection_exists = true;
    } catch(lcio::DataNotAvailableException & e) {
        tot_collection = new LCCollectionVec(lcio::LCIO::TRACKERDATA);
        tot_collection_exists = false;
    }


    // a cell id identifies from which detector component a specific data
    // collection originates.
    // here it is used to encode the telescope / dut plane and the type of
    // data, e.g. zero-suppressed or raw, that is stored in the collection
    CellIDEncoder<TrackerDataImpl>
        cell_encoder(EUTELESCOPE::ZSDATADEFAULTENCODING, mupix_collection);
    if (MUPIX_SENSOR_ID == 601) {
      cell_encoder["sensorID"] = 61;
    }
    else if (MUPIX_SENSOR_ID == 701) {
      cell_encoder["sensorID"] = 71;
    }
    else {
      cell_encoder["sensorID"] = MUPIX_SENSOR_ID;
    }
    cell_encoder["sparsePixelType"] = eutelescope::kEUTelMuPixel;

    CellIDEncoder<TrackerDataImpl>
      cell_encoder_trigger(EUTELESCOPE::ZSDATADEFAULTENCODING, trigger_collection);
    cell_encoder_trigger["sensorID"] = 1;

    CellIDEncoder<TrackerDataImpl>
      cell_encoder_tot(EUTELESCOPE::ZSDATADEFAULTENCODING, tot_collection);
    cell_encoder_tot["sensorID"] = 1;

    // the lcio objects that store the hit and trigger data for a single readout frame
    mupix_frame.reset(new TrackerDataImpl);
    cell_encoder.setCellID(mupix_frame.get());
    trigger_frame.reset(new TrackerDataImpl);
    cell_encoder_trigger.setCellID(trigger_frame.get());
    tot_frame.reset(new TrackerDataImpl);
    cell_encoder_trigger.setCellID(tot_frame.get());

    // a convenience object that encodes the sparse pixel data and external triggers into an
    // eutelescope-specific format and stores them in the given readout frames
    pixels.reset( new EUTelTrackerDataInterfacerImpl<EUTelMuPixel>( mupix_frame.get() ) );
    triggers.reset( new EUTelTrackerDataTriggerInterfacer ( trigger_frame.get() ) );
    tots.reset( new EUTelTrackerDataTriggerInterfacer ( tot_frame.get() ) );


    if (const RawDataEvent * source = dynamic_cast<const RawDataEvent *> (ev)){
      if (const RawDataEvent * next_source = dynamic_cast<const RawDataEvent *> (next)){
        if (const RawDataEvent * nextNext_source = dynamic_cast<const RawDataEvent *> (nextNext)){

	TelescopeFrame data=TelescopeFrame();
	//cout << "current: " << source->NumBlocks() << " next: " << next_source->NumBlocks() << endl;
	for ( uint i = 0; i < source->NumBlocks(); ++i ) {
	  data.from_bytes( source->GetBlock( i ) );

	  /* Extract TLU trigger and other triggers
	     timestamp: RawTrigger::timestamp() returns the absolute timestamp, not the difference
	     between frame and trigger,
	     tag: 0x1 for TLU trigger, 0xBA for "normal" triggers

	     Create ExternalTrigger object which stores timestamps (64 bits) and labels (16 bits)
	     of external triggers,
	     add to triggers LCIO data object
	  */
	  RawTrigger trig;
//	  cout << "size of triggers: " << data.num_triggers() << endl;
	  for ( uint j = 0; j < data.num_triggers(); j++) {
	    trig = data.get_trigger(j);
	    EUTelExternalTrigger eutel_trigger(
					       trig.timestamp(),
					       trig.tag()
					       );
	    triggers->addExternalTrigger( &eutel_trigger );
	  }

	  /* Extract ToT information
	     Create ExternalTrigger object which stores timestamp (48 bits) plus length of ToT (8 bits)
	     and a  label to identify as ToT (0x2),
	     add to triggers LCIO data object
	  */
	  RawTimeOverThreshold tot;
	  for ( uint j = 0; j < data.num_tots(); j++) {
	    tot = data.get_tot(j);
	    EUTelExternalTrigger eutel_tot(
					   (uint64_t)( (tot.timestamp() & 0xFFFFFFFFFFFF)<<8 | (tot.length() & 0xFF ) ),
					   0x2
	  				   );
	    tots->addExternalTrigger( &eutel_tot );
	  }

	  /*
	     Create MuPixel object which also stores the hit and frame timestamp,
	     add to pixels LCIO data object
	  */
	  for ( uint j = 0; j < data.num_hits(); ++j ) {
	    // Caution: columns and rows swapped for EUTelescope analysis
	    // of March 2016 DESY data
	    // since 90 degree rotation is not supported by the framework
	    EUTelMuPixel pixel(
			       //data.get_hit(j,MUPIX_TYPE).column(),
			       //data.get_hit(j,MUPIX_TYPE).row(),
			       data.get_hit(j,MUPIX_TYPE).row(),
			       data.get_hit(j,MUPIX_TYPE).column(),
			       1, // binary
			       0,
			       data.get_hit(j,MUPIX_TYPE).timestamp_raw(),  // hit time stamp (8 bit)
			       data.timestamp() & 0xFFFFFFFF  // frame time stamp, 64 bit!! LCIO can only store floats
			       );
	    if ( data.get_hit(j,MUPIX_TYPE).row() < MUPIX_SENSOR_NUM_COLS && data.get_hit(j,MUPIX_TYPE).column() < MUPIX_SENSOR_NUM_ROWS ) {
	      pixels->addSparsePixel( &pixel );
	    }
	    else
	      cout << "WARNING: col = " << data.get_hit(j,MUPIX_TYPE).row() << ", row = " << data.get_hit(j,MUPIX_TYPE).column() << endl;

	  }
	  data.clear();
	}

	for ( uint i = 0; i < next_source->NumBlocks(); ++i ) {
	  data.from_bytes( next_source->GetBlock( i ) );

	  /* Extract TLU trigger and other triggers
	     timestamp: RawTrigger::timestamp() returns the absolute timestamp, not the difference
	     between frame and trigger,
	     tag: 0x1 for TLU trigger, 0xBA for "normal" triggers

	     Create ExternalTrigger object which stores timestamps (64 bits) and labels (16 bits)
	     of external triggers,
	     add to triggers LCIO data object
	  */
	  RawTrigger trig;
	  for ( uint j = 0; j < data.num_triggers(); j++) {
	    trig = data.get_trigger(j);
	    EUTelExternalTrigger eutel_trigger(
					       trig.timestamp(),
					       trig.tag()
					       );
	    triggers->addExternalTrigger( &eutel_trigger );
	  }

	  /* Extract ToT information
	     Create ExternalTrigger object which stores timestamp (48 bits) plus length of ToT (8 bits)
	     and a  label to identify as ToT (0x2),
	     add to triggers LCIO data object
	  */
	  RawTimeOverThreshold tot;
	  for ( uint j = 0; j < data.num_tots(); j++) {
	    tot = data.get_tot(j);
	    EUTelExternalTrigger eutel_tot(
	  				   (uint64_t)( (tot.timestamp() & 0xFFFFFFFFFFFF)<<8 | (tot.length() & 0xFF ) ),
					   0x2
	  				   );
	    tots->addExternalTrigger( &eutel_tot );
	  }

	  /*
	     Create MuPixel object which also stores the hit and frame timestamp,
	     add to pixels LCIO data object
	  */
	  for ( uint j = 0; j < data.num_hits(); ++j ) {
	    // Caution: columns and rows swapped for EUTelescope analysis
	    // since 90 degree rotation is not supported by the framework
	    EUTelMuPixel pixel(
			       //data.get_hit(j,MUPIX_TYPE).column(),
			       //data.get_hit(j,MUPIX_TYPE).row(),
			       data.get_hit(j,MUPIX_TYPE).row(),
			       data.get_hit(j,MUPIX_TYPE).column(),
			       1, // binary
			       0,
			       data.get_hit(j,MUPIX_TYPE).timestamp_raw(),  // hit time stamp (8 bit)
			       data.timestamp() & 0xFFFFFFFF  // frame time stamp, 64 bit!! LCIO can only store floats
			       );
	    if ( data.get_hit(j,MUPIX_TYPE).row() < MUPIX_SENSOR_NUM_COLS && data.get_hit(j,MUPIX_TYPE).column() < MUPIX_SENSOR_NUM_ROWS ) {
	      pixels->addSparsePixel( &pixel );
	    }
	    else
	      cout << "WARNING: col = " << data.get_hit(j,MUPIX_TYPE).row() << ", row = " << data.get_hit(j,MUPIX_TYPE).column() << endl;

	  }
	  data.clear();
	}

//	third cycle
	for ( uint i = 0; i < nextNext_source->NumBlocks(); ++i ) {
	  data.from_bytes( nextNext_source->GetBlock( i ) );

	  /* Extract TLU trigger and other triggers
	     timestamp: RawTrigger::timestamp() returns the absolute timestamp, not the difference
	     between frame and trigger,
	     tag: 0x1 for TLU trigger, 0xBA for "normal" triggers

	     Create ExternalTrigger object which stores timestamps (64 bits) and labels (16 bits)
	     of external triggers,
	     add to triggers LCIO data object
	  */
	  RawTrigger trig;
	  for ( uint j = 0; j < data.num_triggers(); j++) {
	    trig = data.get_trigger(j);
	    EUTelExternalTrigger eutel_trigger(
					       trig.timestamp(),
					       trig.tag()
					       );
	    triggers->addExternalTrigger( &eutel_trigger );
	  }

	  /* Extract ToT information
	     Create ExternalTrigger object which stores timestamp (48 bits) plus length of ToT (8 bits)
	     and a  label to identify as ToT (0x2),
	     add to triggers LCIO data object
	  */
	  RawTimeOverThreshold tot;
	  for ( uint j = 0; j < data.num_tots(); j++) {
	    tot = data.get_tot(j);
	    EUTelExternalTrigger eutel_tot(
	  				   (uint64_t)( (tot.timestamp() & 0xFFFFFFFFFFFF)<<8 | (tot.length() & 0xFF ) ),
					   0x2
	  				   );
	    tots->addExternalTrigger( &eutel_tot );
	  }

	  /*
	     Create MuPixel object which also stores the hit and frame timestamp,
	     add to pixels LCIO data object
	  */
	  for ( uint j = 0; j < data.num_hits(); ++j ) {
	    // Caution: columns and rows swapped for EUTelescope analysis
	    // since 90 degree rotation is not supported by the framework
	    EUTelMuPixel pixel(
			       //data.get_hit(j,MUPIX_TYPE).column(),
			       //data.get_hit(j,MUPIX_TYPE).row(),
			       data.get_hit(j,MUPIX_TYPE).row(),
			       data.get_hit(j,MUPIX_TYPE).column(),
			       1, // binary
			       0,
			       data.get_hit(j,MUPIX_TYPE).timestamp_raw(),  // hit time stamp (8 bit)
			       data.timestamp() & 0xFFFFFFFF  // frame time stamp, 64 bit!! LCIO can only store floats
			       );
	    if ( data.get_hit(j,MUPIX_TYPE).row() < MUPIX_SENSOR_NUM_COLS && data.get_hit(j,MUPIX_TYPE).column() < MUPIX_SENSOR_NUM_ROWS ) {
	      pixels->addSparsePixel( &pixel );
	    }
	    else
	      cout << "WARNING: col = " << data.get_hit(j,MUPIX_TYPE).row() << ", row = " << data.get_hit(j,MUPIX_TYPE).column() << endl;

	  }
	  data.clear();
	}

	// hand over ownership over the readout frame to the lcio collection
	mupix_collection->push_back( mupix_frame.release() );
	trigger_collection->push_back( trigger_frame.release() );
	tot_collection->push_back( tot_frame.release() );

	if (!mupix_collection_exists && (mupix_collection->size() != 0)) {
	  dest.addCollection(mupix_collection, MUPIX_COLLECTION_NAME);
	}
	else {
	  cout << "FAILED to convert Mupix event to RawDataEvent" << endl;
	}
	if (!trigger_collection_exists && (trigger_collection->size() != 0)) {
	  dest.addCollection(trigger_collection, TRIGGER_COLLECTION_NAME);
	}
	else {
	  cout << "FAILED to convert trigger event to RawDataEvent" << endl;
	}
	if (!tot_collection_exists && (tot_collection->size() != 0)) {
	  dest.addCollection(tot_collection, TOT_COLLECTION_NAME);
	}
	else {
	  cout << "FAILED to convert tot event to RawDataEvent" << endl;
	}
        } //if case for 3rd cycle
      } //if case for 2nd cycle
    } //if case for 1st cycle

    return true;
}
  
#endif  
  
} // namespace eudaq
