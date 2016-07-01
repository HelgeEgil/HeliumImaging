#ifndef TrackCollection_h
#define TrackCollection_h

#include <vector>

#include <TClonesArray.h>

#include "Classes/Track/Track.h"

class Hits;

using namespace std;

struct trackCluster { 
  int track; 
  int cluster;
};

class Tracks : public TObject {

   private:
      TClonesArray tracks_;
      TClonesArray clustersWithoutTrack_;

   public:

      Tracks() : tracks_("Track", 1000), clustersWithoutTrack_("Cluster", 5000) {}
      Tracks(Int_t nTracks) : tracks_("Track", nTracks), clustersWithoutTrack_("Cluster", nTracks*100) {}
      virtual ~Tracks(); 

      virtual void appendTrack(Track *copyTrack, Int_t startOffset = 0);
      virtual void appendClustersWithoutTrack(TClonesArray *clustersWithoutTrack);

      virtual void clearTracks() { tracks_.Clear("C"); }
      virtual void Clear(Option_t * = "");
      virtual void SetOwner(Bool_t val) { tracks_.SetOwner(val); }

      virtual Int_t GetEntriesFast() { return tracks_.GetEntriesFast(); }
      virtual Int_t GetEntriesFastClustersWithoutTrack() { return clustersWithoutTrack_.GetEntriesFast(); }
      virtual Int_t GetEntries() { return tracks_.GetEntries(); }
      virtual Int_t GetEntriesFast(Int_t i) { return At(i)->GetEntriesFast(); }

      Float_t getSinuosity(Int_t i) { return At(i)->getSinuosity(); }
      Float_t getSlopeAngle(Int_t i) { return At(i)->getSlopeAngle(); }
      Float_t getTrackLengthmm(Int_t i) { return At(i)->getTrackLengthmm(); }
      Float_t getTrackScore(Int_t i) { return At(i)->getTrackScore(); }
      TClonesArray * getClustersWithoutTrack() { return (TClonesArray*) &clustersWithoutTrack_; }
		virtual Bool_t isUsedClustersInTrack(Int_t i) { return At(i)->isUsedClustersInTrack(); }
		virtual Int_t getNumberOfConflictClusters(Int_t i) { return At(i)->getNumberOfConflictClusters(); }
		Int_t getTrackIdxFromFirstLayerEID(Int_t eventID);
		vector<Int_t> * getTracksWithConflictClusters();
		vector<Int_t> * getConflictingTracksFromTrack(Int_t trackIdx);
		vector<Int_t> * getTracksFromCluster(Cluster * cluster);
		Int_t getTrackIdxFromCluster(Cluster * cluster);

      virtual void extrapolateToLayer0();
      virtual void splitSharedClusters();
		virtual void matchWithEventIDs(Hits *eventIDs);
      Int_t getClosestCluster(vector<trackCluster> clusters, Cluster* interpolatedCluster);

      virtual Track* At(Int_t i) { return ((Track*) tracks_.At(i)); }

      virtual void removeTrack(Track *t) { tracks_.Remove((TObject*) t); }
      virtual TObject* removeTrackAt(Int_t i) { return tracks_.RemoveAt(i); }
      void sortTrackByLayer(Int_t track);

      void checkLayerOrientation();
      void doFit();
		
		void removeTracksLeavingDetector(); 
		void removeTrackCollisions();
		void retrogradeTrackImprovement(Clusters * clusters);

   ClassDef(Tracks,2);
};
#endif