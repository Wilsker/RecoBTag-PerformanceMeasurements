#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <math.h>
#include <fstream>

using namespace std;

void SommaQuadratura(void){
  const int numbins = 7;
  double ptrange[8]  = {30, 50, 70, 100, 140, 200, 300, 600};

  //Go in the folder with nominal values
  //Backup the original .csv
  //Remove quotation marks from around SF in the .csv file.
  //Use the following command to get a print out of the central values of the SFs:
  //$>cat kin_calib.csv | grep central | awk '{printf "%.12f\n", $11}'

  //DeepCSV SF's from kin_calib.csv
  double sfvals[21]  = {
    0.964403450489,
    0.958738923073,
    0.958408653736,
    0.959267616272,
    0.949074387550,
    0.932909905910,
    0.906148135662,
    0.934341728687,
    0.931151211262,
    0.929023683071,
    0.927854001522,
    0.921115934849,
    0.892432928085,
    0.802657723427,
    0.915498316288,
    0.910543084145,
    0.898578941822,
    0.901119112968,
    0.898955643177,
    0.835006892681,
    0.740537703037
  };


  //Total down/up systematic DeepCSV SF value taken from kin_calib.csv as symmettric anyway.
  double sferrs_up[21]    = {
    0.9459740,
    0.9476040,
    0.9521060,
    0.9503020,
    0.9401600,
    0.9137290,
    0.8385640,
    0.9224420,
    0.9223520,
    0.9227710,
    0.9186130,
    0.9082550,
    0.8640880,
    0.7117250,
    0.9031930,
    0.9002500,
    0.8894290,
    0.8906010,
    0.8818870,
    0.7927960,
    0.6190180
  };

  double sferrs_down[21]    = {
    0.9459740,
    0.9476040,
    0.9521060,
    0.9503020,
    0.9401600,
    0.9137290,
    0.8385640,
    0.9224420,
    0.9223520,
    0.9227710,
    0.9186130,
    0.9082550,
    0.8640880,
    0.7117250,
    0.9031930,
    0.9002500,
    0.8894290,
    0.8906010,
    0.8818870,
    0.7927960,
    0.6190180
  };

  double hdampup[21] = {
    0.977384507656,
    0.965317666531,
    0.964261412621,
    0.964901208878,
    0.956097483635,
    0.956137657166,
    0.944572389126,
    0.943479180336,
    0.932665109634,
    0.930302739143,
    0.932187914848,
    0.921059787273,
    0.898424267769,
    0.871245741844,
    0.924555718899,
    0.915210366249,
    0.895260810852,
    0.901662230492,
    0.890772283077,
    0.844231963158,
    0.796834588051
  };

  double hdampdown[21] = {
    0.974202036858,
    0.960251629353,
    0.956173121929,
    0.954126775265,
    0.948036372662,
    0.939924597740,
    0.873578846455,
    0.936306834221,
    0.925228774548,
    0.928488314152,
    0.920737266541,
    0.917659282684,
    0.911872982979,
    0.811187088490,
    0.915489673615,
    0.909133791924,
    0.894016265869,
    0.892475545406,
    0.897331893444,
    0.861937344074,
    0.723823547363
  };

  double tuneCP5up[21] = {
    0.972653031349,
    0.966113150120,
    0.959789633751,
    0.958514273167,
    0.956540346146,
    0.938754200935,
    0.948832333088,
    0.941490352154,
    0.936455368996,
    0.925165534019,
    0.924914062023,
    0.921061456203,
    0.895212352276,
    0.866632401943,
    0.923632800579,
    0.913440227509,
    0.893759608269,
    0.901933610439,
    0.893895089626,
    0.828745543957,
    0.775218427181
  };

  double tuneCP5down[21] = {
    0.972041070461,
    0.968306899071,
    0.961091816425,
    0.959683477879,
    0.953997552395,
    0.926793694496,
    0.935585260391,
    0.934997737408,
    0.937920153141,
    0.928419113159,
    0.932838022709,
    0.922718524933,
    0.873470962048,
    0.856534123421,
    0.920037806034,
    0.919371783733,
    0.892896950245,
    0.904601275921,
    0.902914583683,
    0.802478313446,
    0.740355134010
  };

  double mtop171p5[21] = {
    0.965444445610,
    0.965464591980,
    0.956845402718,
    0.967020153999,
    0.961651325226,
    0.951460540295,
    0.936724483967,
    0.928632915020,
    0.935040831566,
    0.925760984421,
    0.929636180401,
    0.933834910393,
    0.932835757732,
    0.805638849735,
    0.911275506020,
    0.918061912060,
    0.890947699547,
    0.905605494976,
    0.915576040745,
    0.880870044231,
    0.763223767281
  };

  double mtop173p5[21] = {
    0.970236480236,
    0.960447788239,
    0.956878066063,
    0.955044090748,
    0.948556125164,
    0.929343998432,
    0.936445474625,
    0.935182034969,
    0.931875944138,
    0.925464212894,
    0.920246303082,
    0.912720680237,
    0.884771943092,
    0.819128155708,
    0.916359961033,
    0.915293335915,
    0.897480309010,
    0.888003230095,
    0.890718638897,
    0.824286997318,
    0.799155771732
  };

  double nonttXsecUp[21] = {
    0.972400903702,
    0.963243961334,
    0.961661159992,
    0.961733996868,
    0.951136887074,
    0.934233486652,
    0.899106979370,
    0.942772090435,
    0.935796439648,
    0.932412385941,
    0.930464029312,
    0.922831773758,
    0.892579138279,
    0.790011823177,
    0.923797011375,
    0.915125846863,
    0.901773512363,
    0.903671443462,
    0.900027751923,
    0.833761394024,
    0.725949585438
  };

double nonttXsecDown[21] = {
  0.956089735031,
  0.954091787338,
  0.955052196980,
  0.956717252731,
  0.946930944920,
  0.931380748749,
  0.912804126740,
  0.925609230995,
  0.926365256310,
  0.925532102585,
  0.925156652927,
  0.919348061085,
  0.892279624939,
  0.816439330578,
  0.907074332237,
  0.905797004700,
  0.895254254341,
  0.898445963860,
  0.897880494595,
  0.836398184299,
  0.757309019566
};

  cout<<"DeepCSV discriminator"<<endl;
  for(int wp=0; wp<3; wp++){
    for(int i=0;i<numbins;i++){
      int pos = i+wp*numbins;

      //double errini         = fabs(sfvals[pos]-sferrs[pos]);
      double errini_up         = fabs(sfvals[pos]-sferrs_up[pos]);
      double errini_down         = fabs(sfvals[pos]-sferrs_down[pos]);

      //cout << "sfvals[pos]: " << sfvals[pos] << endl;
      double hdamp_up = hdampup[pos] - sfvals[pos];
      double hdamp_down = hdampdown[pos] - sfvals[pos];
      /*cout << "hdampup[pos]: " << hdampup[pos] << endl;
      cout << "hdampdown[pos]: " << hdampdown[pos] << endl;
      cout << "hdamp_up err: " << hdamp_up << endl;
      cout << "hdamp_down err: " << hdamp_down << endl;*/

      double tuneCP5_up = tuneCP5up[pos] - sfvals[pos];
      double tuneCP5_down = tuneCP5down[pos] - sfvals[pos];
      /*cout << "tuneCP5up[pos]: " << tuneCP5up[pos] << endl;
      cout << "tuneCP5down[pos]: " << tuneCP5down[pos] << endl;
      cout << "tuneCP5_up err: " << tuneCP5_up << endl;
      cout << "tuneCP5_down err: " << tuneCP5_down << endl;*/

      double mtop_175p5 = mtop175p5[pos] - sfvals[pos];
      double mtop_169p5 = mtop169p5[pos] - sfvals[pos];
      /*cout << "mtop175p5[pos]: " << mtop175p5[pos] << endl;
      cout << "mtop169p5[pos]: " << mtop169p5[pos] << endl;
      cout << "mtop_175p5 err: " << mtop_175p5 << endl;
      cout << "mtop_169p5 err: " << mtop_169p5 << endl;*/

      double nonttXsec_Up = nonttXsecUp[pos] - sfvals[pos];
      double nonttXsec_Down = nonttXsecDown[pos] - sfvals[pos];
      /*cout << "nonttXsecUp[pos]: " << nonttXsecUp[pos] << endl;
      cout << "nonttXsecDown[pos]: " << nonttXsecDown[pos] << endl;
      cout << "nonttXsec_Up err: " << nonttXsec_Up << endl;
      cout << "nonttXsec_Down err: " << nonttXsec_Down << endl;*/


      //New total error and values
      //double toterr = sqrt( pow(errini,2) + pow(hdamp_err,2) + pow(tuneCP5_err,2) + pow(mtop_err,2) + pow(nonttXsec_err,2) );
      double toterr_up = sqrt( pow(errini_up,2) + pow(hdamp_up,2) + pow(tuneCP5_up,2) + pow(mtop_175p5,2) + pow(nonttXsec_Up,2) );
      double toterr_down = sqrt( pow(errini_down,2) + pow(hdamp_down,2) + pow(tuneCP5_down,2) + pow(mtop_169p5,2) + pow(nonttXsec_Down,2) );
      double valueup = sfvals[pos]+toterr_up;
      double valuedown = sfvals[pos]-toterr_down;

      //cout << "wp: " << wp << ", jet pt range: " << ptrange[i] << " - " << ptrange[i+1] << " , toterr: " << toterr << endl;
      //cout << "wp: " << wp << ", jet pt range: " << ptrange[i] << " - " << ptrange[i+1] << " , hdamp_err: " << hdamp_err << endl;
      //cout << "wp: " << wp << ", jet pt range: " << ptrange[i] << " - " << ptrange[i+1] << " , tuneCP5_err: " << tuneCP5_err << endl;
      //cout << "wp: " << wp << ", jet pt range: " << ptrange[i] << " - " << ptrange[i+1] << " , mtop_err: " << mtop_err << endl;
      //cout << "wp: " << wp << ", jet pt range: " << ptrange[i] << " - " << ptrange[i+1] << " , nonttXsec_err: " << nonttXsec_err << endl;


      //cout<<wp<<", kin, down_hdamp, 0, -2.4, 2.4, "<<ptrange[i]<<", "<<ptrange[i+1]<<", 0, 1, \""<<hdampdown[pos]<<"\""<<endl;
      //cout<<wp<<", kin, up_hdamp, 0, -2.4, 2.4, "<<ptrange[i]<<", "<<ptrange[i+1]<<", 0, 1, \""<<hdampup[pos]<<"\""<<endl;
      //cout<<wp<<", kin, down_tuneCP5, 0, -2.4, 2.4, "<<ptrange[i]<<", "<<ptrange[i+1]<<", 0, 1, \""<<tuneCP5down[pos] <<"\""<<endl;
      //cout<<wp<<", kin, up_tuneCP5, 0, -2.4, 2.4, "<<ptrange[i]<<", "<<ptrange[i+1]<<", 0, 1, \""<< tuneCP5up[pos] <<"\""<<endl;
      //cout<<wp<<", kin, down_mtop, 0, -2.4, 2.4, "<<ptrange[i]<<", "<<ptrange[i+1]<<", 0, 1, \""<<mtop169p5[pos]<<"\""<<endl;
      //cout<<wp<<", kin, up_mtop, 0, -2.4, 2.4, "<<ptrange[i]<<", "<<ptrange[i+1]<<", 0, 1, \""<< mtop175p5[pos] <<"\""<<endl;
      //cout<<wp<<", kin, up_nonttXsec, 0, -2.4, 2.4, "<<ptrange[i]<<", "<<ptrange[i+1]<<", 0, 1, \""<<nonttXsecUp[pos]<<"\""<<endl;
      //cout<<wp<<", kin, down_nonttXsec, 0, -2.4, 2.4, "<<ptrange[i]<<", "<<ptrange[i+1]<<", 0, 1, \""<<nonttXsecDown[pos]<<"\""<<endl;
      cout<<wp<<", kin, down, 0, -2.4, 2.4, "<<ptrange[i]<<", "<<ptrange[i+1]<<", 0, 1, \""<<valuedown<<"\""<<endl;
      //cout<<wp<<", kin, up, 0, -2.4, 2.4, "<<ptrange[i]<<", "<<ptrange[i+1]<<", 0, 1, \""<<valueup<<"\""<<endl;

    }
  }
}
