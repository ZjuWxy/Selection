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

public class entry_2 {

    /**
     *  Characteristic parameters of LC-MS
     */
    double TOTALmz;
    double[] eachMZ;
    double[] eachTIC;
    double[] eachRT;

    /**
     *  Constructor function
     */
    public entry_2(double TOTALmz_IN, double[] eachMZ_IN, double[] eachTIC_IN, double[] eachRT_IN) {
        TOTALmz = TOTALmz_IN;
        eachMZ = eachMZ_IN;
        eachTIC = eachTIC_IN;
        eachRT = eachRT_IN;
    }
}











