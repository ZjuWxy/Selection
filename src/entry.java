/**
 *  This class provides related parameters of LC-MS data
 *  and counting the #number of the possible products.
 *  <p>
 *  The halogenated products were screened from many aspects.
 *  <p>
 *  For additional documentation, see
 *  https://github.com/ZjuWxy/SelectionForDBPs
 *  @author Wang Xinyang
 */

public class entry {

    /**
     *  Characteristic parameters of LC-MS
     */
    int index;            // sample subscript
    double RT;            // retention time
    double mz;            // mass-to-charge ratio
    double MassDefect;    // quality defect
    double TIC;           // peak intensity
    int NumMerged;
    int Quality;

    /**
     *  Constructor function
     */
    public entry(int indexIN, double RTIN, double mzIN, double MassDefectIN,
         double TICIN, int NumMergedIN, int QualityIN) {
        index = indexIN;            // read sample subscript
        RT = RTIN;                  // read retention time
        mz = mzIN;                  // read mass-to-charge ratio
        MassDefect = MassDefectIN;  // read quality defect
        TIC = TICIN;                // read peak intensity
        NumMerged = NumMergedIN;
        Quality = QualityIN;
    }
}











