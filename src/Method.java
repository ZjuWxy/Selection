import java.io.FileNotFoundException;
import java.io.PrintStream;
import static java.lang.Math.abs;

/**
 *  This class provides methods for screening of halogenated nucleotides
 *  and counting the #number of the possible products.
 *  <p>
 *  I sort the samples by mass charge ratio from smallest to largest with PeakView software.
 *  The identification of unknown halogenated products requires consideration of many factors,
 *  such as mass spectrum signal intensity, mass/charge ratio, cracking mode, etc.,
 *  which are crucial for the detection and identification of halogenated products.
 *  <p>
 *  According to these characteristics, write the relevant code algorithm,
 *  using the Java programming language,
 *  complete the input and output of the data file,
 *  and the above various factors and the resolution and accuracy of the mass spectrometer.
 *  <p>
 *  For additional documentation, see
 *  https://github.com/ZjuWxy/Selection
 *
 *  @author Wang Xinyang
 */

public class Method {
    /**
     * Create a list of input fileName
     * and the relevant data length.
     */
    static String[] inputName = {"AMP-NEG.txt", "CMP-NEG.txt", "UMP-NEG.txt", "GMP-NEG.txt"};
    static int[] len = new int[inputName.length];
    static int[][] count = new int[inputName.length][5];
    static int[][] count_78_96 = new int[inputName.length][5];

    /**
     * parameters M/Z,
     */
    private static double mzDiff = 2;

    /**
     * parameters TIC corresponding 1-Cl, 2-Cl, 1-Br, 2-Br, 1-Br-1-Cl
     */
    static double[] TICDiff = {32.0, 63.9, 97.4, 194.6, 97.3}; // %

    /**
     * Decide relevant parameter tolerance,
     * including M/Z, TIC and RT.
     */
    private static double limMZ  = 0.01;
    private static double limTIC = 12; // %
    private static double limRT  = 0.15;


    /**
     * Read in date from PeakView list.
     */
    public static entry_1[] readFromTxt (String st, int inputI){
        // create new objects
        In in = new In(st);
        len[inputI] = in.readInt();
        entry_1[] res = new entry_1[len[inputI]+1];

        // read in all feature
        for (int i = 1; i <= len[inputI]; i++)
            res[i] = new entry_1(in.readInt(),in.readDouble(),in.readDouble(),
                    in.readDouble(),in.readDouble(), in.readInt(), in.readInt());
        return res;
    }

    static void specialPrint(entry_1[] entries, String outputName, int TICI, int inputI, double minMZ) throws FileNotFoundException {
        PrintStream ps = new PrintStream(outputName);;
        System.setOut(ps);

        // lLim & rLim
        double lMZ = mzDiff - limMZ;
        double rMZ = mzDiff + limMZ;
        double lTIC = TICDiff[TICI] - limTIC;
        double rTIC = TICDiff[TICI] + limTIC;

        // print input & output file name
        System.out.println("Input  file name: " + inputName[inputI]);
        System.out.println("Output file name: " + outputName);
        System.out.println();

        /**
         * Decide relevant parameter tolerance,
         * including M/Z, TIC and RT.
         */
        System.out.println("1.参数设置: ");
        System.out.println("   样本数量: " + len[inputI]);
        System.out.println("   m/z偏差: " + limMZ);
        System.out.println("   TIC偏差: " + limTIC + "%");
        System.out.println("   RT偏差 : " + limRT);
        System.out.println();
        System.out.println("2.ChemDraw标准: ");
        System.out.println("   m/z差: " + mzDiff);
        System.out.println("   TIC: 100.0% , " + String.format("%.1f", TICDiff[TICI]) + "%");
        System.out.println();

        System.out.println("3.筛选标准: ");
        System.out.println("   质荷比:  M : (M+" + lMZ + ") ~ (M+" + rMZ + ") ");
        System.out.println("   丰度比:  100.0% , " + String.format("%.1f", lTIC) + "%" + " ~ "
                + String.format("%.1f", rTIC) + "%");
        System.out.println("   保留时间差:  |RT_1 - RT_2| < " + limRT);
        System.out.println();


        // #sample
        count[inputI][TICI] = 0;
        int[] ans_1 = new int[len[inputI] + 1];
        int[] ans_2 = new int[len[inputI] + 1];
        int[] ans_3 = new int[len[inputI] + 1];
        // two iteration
        for (int i = 1; i <= len[inputI]; i++) {
            if(entries[i].mz < minMZ) continue;
            // relevant standard
            double st1, st2, st3;
            double tmp = 100000;
            for (int k = i; k < entries.length; k++) {

                st1 = abs(entries[k].mz - entries[i].mz - mzDiff);
                st2 = abs(entries[k].TIC / entries[i].TIC * 100 - TICDiff[TICI]);
                st3 = abs(entries[i].RT - entries[k].RT);

                double st4, st5, st6;
                for(int j = i; j < entries.length; j++){
                    st4 = abs(entries[j].mz - entries[i].mz - mzDiff);
                    st5 = abs(entries[j].TIC / entries[i].TIC * 100 - 32.0);
                    st6 = abs(entries[j].RT - entries[i].RT);

                    // standards to meet
                    if (st1 <= limMZ && st2 <= limTIC && st3 <= limRT
                            && st4 <= limMZ && st5 <= limTIC && st6 <= limRT) {
                        if(tmp == 100000){
                            count[inputI][TICI] += 1; // 样本计数
                        }
                        if(st2 < tmp){
                            tmp = st2;
                            ans_1[count[inputI][TICI]] = i;
                            ans_2[count[inputI][TICI]] = k;
                            ans_3[count[inputI][TICI]] = j;
                        }
                    }
                }
            }
        }

        // print #sample
        System.out.println("4.可能的卤代样本: ");
        System.out.println("   样本总数: " + len[inputI]);
        System.out.println("   可能的卤代样本总数: " + count[inputI][TICI]);
        System.out.println();

        // creat match object
        Match myTest = null;
        entry_2[] myEntry2 = myTest.readFromMgf(inputName[inputI].substring(0, 3) + "-1-NEG.mgf", len[inputI]);


        // print details
        for (int i = 1; i <= count[inputI][TICI]; i++) {
            System.out.println("   样本序号: " + i);
            System.out.println("   索引号: " + "#" + entries[ans_1[i]].index + " 和 "
                    + "#" + entries[ans_2[i]].index);

            // all features
            System.out.print("   质荷比: " + String.format("%.4f 和 ", entries[ans_1[i]].mz) + String.format("%.4f", entries[ans_2[i]].mz));
            System.out.println();
            System.out.print("   丰度比: 100.0% 和 " + String.format("%.1f", entries[ans_2[i]].TIC / entries[ans_1[i]].TIC * 100.0) + "%");
            System.out.println();
            System.out.print("   保留时间: " + String.format("%.4f 和 ", entries[ans_1[i]].RT) + String.format("%.4f", entries[ans_2[i]].RT));
            System.out.println();

            /**
             * The output of each selected mass-charge ratio peak.
             */
            System.out.println("   #" + entries[ans_1[i]].index +
                    "  质荷比: " + String.format("%.4f", entries[ans_1[i]].mz) +
                    ";   峰强度: " + entries[ans_1[i]].TIC +
                    ";   保留时间: " + String.format("%.4f", entries[ans_1[i]].RT));
            System.out.println("   #" + entries[ans_2[i]].index +
                    "  质荷比: " + String.format("%.4f", entries[ans_2[i]].mz) +
                    ";   峰强度: " + entries[ans_2[i]].TIC +
                    ";   保留时间: " + String.format("%.4f", entries[ans_2[i]].RT));
            System.out.println("   #" + entries[ans_3[i]].index +
                    "  质荷比: " + String.format("%.4f", entries[ans_3[i]].mz) +
                    ";   峰强度: " + entries[ans_3[i]].TIC +
                    ";   保留时间: " + String.format("%.4f", entries[ans_3[i]].RT));

            // MS/MS
            int flag_78_1 = 0, flag_96_1 = 0;
            int flag_78_2 = 0, flag_96_2 = 0;

            System.out.println("   #" + entries[ans_1[i]].index + "的碎片离子：");

            for(int j = 1; j <= len[inputI]-10; j++){
                double diff = abs(myEntry2[j].TOTALmz - entries[ans_1[i]].mz);
                int inner_96_1 = 0, inner_78_1 = 0;
                if(diff < 0.01){
                    int k_pre = -1;
                    for(int k = 0; k <= 400; k++){
                        if(myEntry2[j].eachMZ[k] == 0) break;
                        double mz_Diff_78 = abs(myEntry2[j].eachMZ[k] - 78.96);
                        double mz_Diff_96 = abs(myEntry2[j].eachMZ[k] - 96.97);
                        if(mz_Diff_78 < 0.01){
                            inner_78_1 = 1;
                            k_pre = k;
                        }
                        if(mz_Diff_96 < 0.01){
                            inner_96_1 = 1;
                        }

                        if(inner_78_1 == 1 && inner_96_1 == 1){
                            System.out.println("   " + myEntry2[j].eachMZ[k_pre]);
                            System.out.println("   " + myEntry2[j].eachMZ[k]);
                            flag_78_1 = 1;
                            flag_96_1 = 1;
                            break;
                        }
                    }
                }
                if(flag_78_1 == 1 && flag_96_1 == 1) break;
            }
            if(flag_78_1 == 0 && flag_96_1 == 0){
                System.out.println("   无特征碎片离子");
            }



            System.out.println("   #" + entries[ans_2[i]].index + "的碎片离子：");

            for(int j = 1; j <= len[inputI]-10; j++){
                double diff = abs(myEntry2[j].TOTALmz - entries[ans_2[i]].mz);
                int inner_96_2 = 0, inner_78_2 = 0;
                if(diff < 0.01){
                    int k_pre = -1;
                    for(int k = 0; k <= 400; k++){
                        if(myEntry2[j].eachMZ[k] == 0) break;
                        double mz_Diff_78 = abs(myEntry2[j].eachMZ[k] - 78.96);
                        double mz_Diff_96 = abs(myEntry2[j].eachMZ[k] - 96.97);
                        if(mz_Diff_78 < 0.01){
                            inner_78_2 = 1;
                            k_pre = k;
                        }
                        if(mz_Diff_96 < 0.01){
                            inner_96_2 = 1;
                        }
                        if(inner_78_2 == 1 && inner_96_2 == 1){
                            System.out.println("   " + myEntry2[j].eachMZ[k_pre]);
                            System.out.println("   " + myEntry2[j].eachMZ[k]);
                            flag_78_2 = 1;
                            flag_96_2 = 1;
                            break;
                        }
                    }
                    if(flag_78_2 == 1 && flag_96_2 == 1) break;
                }
            }
            if(flag_78_2 == 0 && flag_96_2 == 0){
                System.out.println("   无特征碎片离子");
            }
            if(flag_78_1 == 1 && flag_96_1 == 1 && flag_78_2 == 1 && flag_96_2 == 1){
                count_78_96[inputI][TICI] += 1;
            }

            System.out.println();
        }
        System.out.println("   包含碎片78.96 & 96.97 总数：" + count_78_96[inputI][TICI]);
    }

    /**
     * my print function for print stream
     */
    static void myPrint(entry_1[] entries, String outputName, int TICI, int inputI, double minMZ) throws FileNotFoundException {
        if(TICI == 4){
            specialPrint(entries, outputName, TICI, inputI, minMZ);
            return;
        }

        PrintStream ps = new PrintStream(outputName);
        System.setOut(ps);

        // lLim & rLim
        double lMZ = mzDiff - limMZ;
        double rMZ = mzDiff + limMZ;
        double lTIC = TICDiff[TICI] - limTIC;
        double rTIC = TICDiff[TICI] + limTIC;


        // print input & output file name
        System.out.println("Input  file name: " + inputName[inputI]);
        System.out.println("Output file name: " + outputName);
        System.out.println();

        /**
         * Decide relevant parameter tolerance,
         * including M/Z, TIC and RT.
         */
        System.out.println("1.参数设置: ");
        System.out.println("   样本数量: " + len[inputI]);
        System.out.println("   m/z偏差: " + limMZ);
        System.out.println("   TIC偏差: " + limTIC + "%");
        System.out.println("   RT偏差 : " + limRT);
        System.out.println();
        System.out.println("2.ChemDraw标准: ");
        System.out.println("   m/z差: " + mzDiff);
        System.out.println("   TIC: 100.0% , " + String.format("%.1f", TICDiff[TICI]) + "%");
        System.out.println();

        System.out.println("3.筛选标准: ");
        System.out.println("   质荷比:  M : (M+" + lMZ + ") ~ (M+" + rMZ + ") ");
        System.out.println("   丰度比:  100.0% , " + String.format("%.1f", lTIC) + "%" + " ~ "
                + String.format("%.1f", rTIC) + "%");
        System.out.println("   保留时间差:  |RT_1 - RT_2| < " + limRT);
        System.out.println();


        // #sample
        count[inputI][TICI] = 0;
        int[] ans_1 = new int[len[inputI] + 1];
        int[] ans_2 = new int[len[inputI] + 1];
        // two iteration
        for (int i = 1; i <= len[inputI]; i++) {
            if(entries[i].mz < minMZ) continue;
            if(i > 1 && ans_2[count[inputI][TICI]] != 0 && abs(entries[i].mz - entries[ans_2[count[inputI][TICI]]].mz) < 1) continue;
            // relevant standard
            double st1, st2, st3;
            double tmp = 100000;
            for (int k = 1; k < i; k++) {

                st1 = abs(entries[i].mz - entries[k].mz - mzDiff);
                st2 = abs(entries[i].TIC / entries[k].TIC * 100 - TICDiff[TICI]);
                st3 = abs(entries[i].RT - entries[k].RT);

                // standards to meet
                if (st1 <= limMZ && st2 <= limTIC && st3 <= limRT) {
                    if(tmp == 100000){
                        count[inputI][TICI] += 1; // 样本计数
                    }
                    if(st2 < tmp){
                        tmp = st2;
                        ans_1[count[inputI][TICI]] = k;
                        ans_2[count[inputI][TICI]] = i;
                    }
                }
            }
        }

        // print #sample
        System.out.println("4.可能的卤代样本: ");
        System.out.println("   样本总数: " + len[inputI]);
        System.out.println("   可能的卤代样本总数: " + count[inputI][TICI]);
        System.out.println();

        // creat match object
        Match myTest = null;
        entry_2[] myEntry2 = myTest.readFromMgf(inputName[inputI].substring(0, 3) + "-1-NEG.mgf", len[inputI]);

        // print details
        for (int i = 1; i <= count[inputI][TICI]; i++) {
            System.out.println("   样本序号: " + i);
            System.out.println("   索引号: " + "#" + entries[ans_1[i]].index + " 和 "
                    + "#" + entries[ans_2[i]].index);

            // all features
            System.out.print("   质荷比: " + String.format("%.4f 和 ", entries[ans_1[i]].mz) + String.format("%.4f", entries[ans_2[i]].mz));
            System.out.println();
            System.out.print("   丰度比: 100.0% 和 " + String.format("%.1f", entries[ans_2[i]].TIC / entries[ans_1[i]].TIC * 100.0) + "%");
            System.out.println();
            System.out.print("   保留时间: " + String.format("%.4f 和 ", entries[ans_1[i]].RT) + String.format("%.4f", entries[ans_2[i]].RT));
            System.out.println();

            /**
             * The output of each selected mass-charge ratio peak.
             */
            System.out.println("   #" + entries[ans_1[i]].index +
                    "  质荷比: " + String.format("%.4f", entries[ans_1[i]].mz) +
                    ";   峰强度: " + entries[ans_1[i]].TIC +
                    ";   保留时间: " + String.format("%.4f", entries[ans_1[i]].RT));
            System.out.println("   #" + entries[ans_2[i]].index +
                    "  质荷比: " + String.format("%.4f", entries[ans_2[i]].mz) +
                    ";   峰强度: " + entries[ans_2[i]].TIC +
                    ";   保留时间: " + String.format("%.4f", entries[ans_2[i]].RT));

            // MS/MS
            int flag_78_1 = 0, flag_96_1 = 0;
            int flag_78_2 = 0, flag_96_2 = 0;

            System.out.println("   #" + entries[ans_1[i]].index + "的碎片离子：");

            for(int j = 1; j <= len[inputI]-10; j++){
                double diff = abs(myEntry2[j].TOTALmz - entries[ans_1[i]].mz);
                int inner_96_1 = 0, inner_78_1 = 0;
                if(diff < 0.01){
                    int k_pre = -1;
                    for(int k = 0; k <= 400; k++){
                        if(myEntry2[j].eachMZ[k] == 0) break;
                        double mz_Diff_78 = abs(myEntry2[j].eachMZ[k] - 78.96);
                        double mz_Diff_96 = abs(myEntry2[j].eachMZ[k] - 96.97);
                        if(mz_Diff_78 < 0.01){
                            inner_78_1 = 1;
                            k_pre = k;
                        }
                        if(mz_Diff_96 < 0.01){
                            inner_96_1 = 1;
                        }

                        if(inner_78_1 == 1 && inner_96_1 == 1){
                            System.out.println("   " + myEntry2[j].eachMZ[k_pre]);
                            System.out.println("   " + myEntry2[j].eachMZ[k]);
                            flag_78_1 = 1;
                            flag_96_1 = 1;
                            break;
                        }
                    }
                }
                if(flag_78_1 == 1 && flag_96_1 == 1) break;
            }
            if(flag_78_1 == 0 && flag_96_1 == 0){
                System.out.println("   无特征碎片离子");
            }



            System.out.println("   #" + entries[ans_2[i]].index + "的碎片离子：");

            for(int j = 1; j <= len[inputI]-10; j++){
                double diff = abs(myEntry2[j].TOTALmz - entries[ans_2[i]].mz);
                int inner_96_2 = 0, inner_78_2 = 0;
                if(diff < 0.01){
                    int k_pre = -1;
                    for(int k = 0; k <= 400; k++){
                        if(myEntry2[j].eachMZ[k] == 0) break;
                        double mz_Diff_78 = abs(myEntry2[j].eachMZ[k] - 78.96);
                        double mz_Diff_96 = abs(myEntry2[j].eachMZ[k] - 96.97);
                        if(mz_Diff_78 < 0.01){
                            inner_78_2 = 1;
                            k_pre = k;
                        }
                        if(mz_Diff_96 < 0.01){
                            inner_96_2 = 1;
                        }
                        if(inner_78_2 == 1 && inner_96_2 == 1){
                            System.out.println("   " + myEntry2[j].eachMZ[k_pre]);
                            System.out.println("   " + myEntry2[j].eachMZ[k]);
                            flag_78_2 = 1;
                            flag_96_2 = 1;
                            break;
                        }
                    }
                    if(flag_78_2 == 1 && flag_96_2 == 1) break;
                }
            }
            if(flag_78_2 == 0 && flag_96_2 == 0){
                System.out.println("   无特征碎片离子");
            }
            if(flag_78_1 == 1 && flag_96_1 == 1 && flag_78_2 == 1 && flag_96_2 == 1){
                count_78_96[inputI][TICI] += 1;
            }


            System.out.println();
        }

        System.out.println("   包含碎片78.96 & 96.97 总数：" + count_78_96[inputI][TICI]);
    }

    /**
     * The final data were counted,
     * including the number of original samples
     * and the number of possible halogen samples selected.
     */
    static void myDataStatistics() throws FileNotFoundException {
        PrintStream ps = new PrintStream("Data_statistics.txt");;
        System.setOut(ps);
        System.out.println("Data_statistics: ");

        for(int i = 0; i < inputName.length; i++){
            int total = count[i][0] + count[i][1] + count[i][2] + count[i][3];
            int total_2_78_96 = count_78_96[i][0] + count_78_96[i][1] + count_78_96[i][2] + count_78_96[i][3];
            if (inputName[i] == "AMP-NEG.txt") {
                /**
                 * The final data were counted,
                 * including the number of original samples
                 * and the number of possible halogen samples selected.
                 */

                System.out.println("1.AMP: ");
                System.out.println("    AMP-1-Cl     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][0]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][0]);
                System.out.println("    AMP-2-Cl     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][1]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][1]);
                System.out.println("    AMP-1-Br     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][2]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][2]);
                System.out.println("    AMP-2-Br     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][3]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][3]);
                System.out.println("    AMP-1-Br-1-Cl: 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][4]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][4]);
                System.out.println("    第一轮总数：" + total);
                System.out.println("    第二轮总数：" + total_2_78_96);
                System.out.println();
            }
            if (inputName[i] == "CMP-NEG.txt") {
                System.out.println("2.CMP: ");
                System.out.println("    CMP-1-Cl     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][0]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][0]);
                System.out.println("    CMP-2-Cl     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][1]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][1]);
                System.out.println("    CMP-1-Br     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][2]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][2]);
                System.out.println("    CMP-2-Br     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][3]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][3]);
                System.out.println("    CMP-1-Br-1-Cl: 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][4]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][4]);
                System.out.println("    第一轮总数：" + total);
                System.out.println("    第二轮总数：" + total_2_78_96);
                System.out.println();
            }
            if (inputName[i] == "UMP-NEG.txt") {
                System.out.println("3.UMP: ");
                System.out.println("    UMP-1-Cl     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][0]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][0]);
                System.out.println("    UMP-2-Cl     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][1]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][1]);
                System.out.println("    UMP-1-Br     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][2]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][2]);
                System.out.println("    UMP-2-Br     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][3]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][3]);
                System.out.println("    UMP-1-Br-1-Cl: 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][4]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][4]);
                System.out.println("    第一轮总数：" + total);
                System.out.println("    第二轮总数：" + total_2_78_96);
                System.out.println();
            }
            if (inputName[i] == "GMP-NEG.txt") {
                System.out.println("4.GMP: ");
                System.out.println("    GMP-1-Cl     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][0]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][0]);
                System.out.println("    GMP-2-Cl     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][1]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][1]);
                System.out.println("    GMP-1-Br     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][2]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][2]);
                System.out.println("    GMP-2-Br     : 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][3]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][3]);
                System.out.println("    GMP-1-Br-1-Cl: 总样本数-" + len[i] + "  " + "可能的取代产物数-" + count[i][4]);
                System.out.println("    包含碎片78.96 & 96.97 总数：" + count_78_96[i][4]);
                System.out.println("    第一轮总数：" + total);
                System.out.println("    第二轮总数：" + total_2_78_96);
                System.out.println();
            }
        }

        // for each halide
        int sum_2_1Cl = 0, sum_2_2Cl = 0;
        int sum_2_1Br = 0, sum_2_2Br = 0;
        int sum_2_1Br1Cl = 0;

        int sum_1_1Cl = 0, sum_1_2Cl = 0;
        int sum_1_1Br = 0, sum_1_2Br = 0;
        int sum_1_1Br1Cl = 0;

        for(int i = 0; i < count.length; i++){
            sum_1_1Cl += count[i][0];
            sum_1_2Cl += count[i][1];
            sum_1_1Br += count[i][2];
            sum_1_2Br += count[i][3];
            sum_1_1Br1Cl += count[i][4];


            sum_2_1Cl += count_78_96[i][0];
            sum_2_2Cl += count_78_96[i][1];
            sum_2_1Br += count_78_96[i][2];
            sum_2_2Br += count_78_96[i][3];
            sum_2_1Br1Cl += count_78_96[i][4];

        }
        int total_1 = sum_1_1Br + sum_1_2Br + sum_1_1Cl + sum_1_2Cl;
        int total_2 = sum_2_1Br + sum_2_2Br + sum_2_1Cl + sum_2_2Cl;
        // total sum

        System.out.println("第一轮：");
        System.out.println("一氯: " + sum_1_1Cl);
        System.out.println("二氯: " + sum_1_2Cl);
        System.out.println("一溴: " + sum_1_1Br);
        System.out.println("二溴: " + sum_1_2Br);
        System.out.println("一溴一氯: " + sum_1_1Br1Cl);
        System.out.println("总数: " + total_1);
        System.out.println();

        System.out.println("第二轮：");
        System.out.println("一氯: " + sum_2_1Cl);
        System.out.println("二氯: " + sum_2_2Cl);
        System.out.println("一溴: " + sum_2_1Br);
        System.out.println("二溴: " + sum_2_2Br);
        System.out.println("一溴一氯: " + sum_2_1Br1Cl);
        System.out.println("总数: " + total_2);
        System.out.println();
    }

    /**
     * The main function includes one iteration
     * which scans the TXT files of AMP, CMP, UMP and GMP
     * and output the 1-Cl, 2-Cl, 1-Br, 2-Br and 1-Br-1-Cl
     * results to the appropriate TXT files.
     */
    public static void main(String[] args) throws FileNotFoundException {
        for (int i = 0; i < inputName.length; i++) {
            entry_1[] entries = readFromTxt(inputName[i], i); // input txt name
            if (inputName[i] == "AMP-NEG.txt") {
                myPrint(entries, "AMP/AMP-1-Cl.txt", 0, i, 347);
                myPrint(entries, "AMP/AMP-2-Cl.txt", 1, i, 347);
                myPrint(entries, "AMP/AMP-1-Br.txt", 2, i, 347);
                myPrint(entries, "AMP/AMP-2-Br.txt", 3, i, 347);
                myPrint(entries, "AMP/AMP-1-Cl-1-Br.txt", 4, i, 347);
            }
            if (inputName[i] == "CMP-NEG.txt") {
                myPrint(entries, "CMP/CMP-1-Cl.txt", 0, i, 323);
                myPrint(entries, "CMP/CMP-2-Cl.txt", 1, i, 323);
                myPrint(entries, "CMP/CMP-1-Br.txt", 2, i, 323);
                myPrint(entries, "CMP/CMP-2-Br.txt", 3, i, 323);
                myPrint(entries, "CMP/CMP-1-Cl-1-Br.txt", 4, i, 323);
            }
            if (inputName[i] == "UMP-NEG.txt") {
                myPrint(entries, "UMP/UMP-1-Cl.txt", 0, i, 324);
                myPrint(entries, "UMP/UMP-2-Cl.txt", 1, i, 324);
                myPrint(entries, "UMP/UMP-1-Br.txt", 2, i, 324);
                myPrint(entries, "UMP/UMP-2-Br.txt", 3, i, 324);
                myPrint(entries, "UMP/UMP-1-Cl-1-Br.txt", 4, i, 324);
            }
            if (inputName[i] == "GMP-NEG.txt") {
                myPrint(entries, "GMP/GMP-1-Cl.txt", 0, i, 363);
                myPrint(entries, "GMP/GMP-2-Cl.txt", 1, i, 363);
                myPrint(entries, "GMP/GMP-1-Br.txt", 2, i, 363);
                myPrint(entries, "GMP/GMP-2-Br.txt", 3, i, 363);
                myPrint(entries, "GMP/GMP-1-Cl-1-Br.txt", 4, i, 363);
            }
        }

        // final data statistics
        myDataStatistics();
    }
}
