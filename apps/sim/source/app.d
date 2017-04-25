module app;

import std.algorithm;
import std.format;
import std.json;
import std.math;
import std.meta;
import std.path;
import std.range;
import std.stdio;

import dffdd.utils.unit;

import models;
import simmain;


void mainJob()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    auto taskList = new MultiTaskList();

    // ADC&IQ&PA
    foreach(methodName; AliasSeq!(/*"IQISICFHF_X", "WLFHF_LS",*/ "SFHF_LS", "S2FHF_LS", /*"OPH_LS", "WLFHF_LS", "WL_LS", "C2DCMFHF_LS", "P2DCMFHF_LS", *//*"SFHF_RLS", "WL_RLS",*/ /*"OPH_LS",*/ /*"SFHF_LS",*//* "C1DCMFHF_LS",*/ /*"C2SDCMFHF_LS",*/ /*"P1DCMFHF_LS", *//*"P2SDCMFHF_LS",*/
                // "FHF_LMS", "FHF_LS", "OPH_LS", "OPH_RLS", "OPH_LMS", "OCH_LS", "OCH_RLS", "OCH_LMS", "WL_LS", "WL_RLS", "WL_LMS", "L_LS", "L_RLS", "L_LMS" /*"FHF", "PH"*//*, "OPH", "OPHDCM", "OCH", "WL", "L",*/ /*"OPHDCM"*/
            ))
        foreach(learningSymbols; iota(60, 65, 5)) foreach(orthTrainingSymbols; [10000])
        {
            Model[] models;
            string[] dirs;

            foreach(inr; iota(20, 85, 5)) foreach(txp; iota(15, 20, 5))
            foreach(gamma; /*iota(0, 8, 2)*/ [6])
            foreach(beta; /*iota(0, 30, 5)*/ [25])
            {
                Model model;
                model.SNR = 11.dB;
                model.INR = inr.dB;
                model.pa.TX_POWER = txp.dBm;

                // model.withSI = false;
                // model.withSIC = false;

                // 再現する非線形性の選択
                model.useDTXIQ = false;
                model.useDTXPN = false;
                model.useDTXPA = false;
                model.useSTXIQ = true;
                model.useSTXPN = false;
                model.useSTXPA = true;
                model.useSTXIQ2 = false;
                model.useSRXLN = true;
                model.useSRXIQ = true;
                model.useSRXQZ = true;

                model.orthogonalizer.numOfTrainingSymbols = orthTrainingSymbols;

                if(methodName[0] == 'O')
                    model.orthogonalizer.enabled = true;
                else
                    model.orthogonalizer.enabled = false;

                // ベースバンド信号波形の出力
                model.outputWaveform = false;

                model.channel.taps = 64;
                model.firFilter.taps = 64;

                model.outputBER = false;

              static if(methodName.canFind("DCM"))
              {
                model.learningSymbols = 10;
                model.learningCount = 3;
              }
              else
              {
                model.learningSymbols = learningSymbols;
                model.learningCount = 1;
              }

              static if(methodName.split("_")[0].endsWith("FHF"))
              {
                model.swappedSymbols = 100000;
                model.numOfFilterTrainingSymbols = 1000;
              }
              else
              {
                model.swappedSymbols = 0;
                model.numOfFilterTrainingSymbols = 1000;
              }

                model.basisFuncsSelection.imageMargin = (-beta).dB;
                model.basisFuncsSelection.noiseMargin = gamma.dB;

                models ~= model;

              static if(methodName.endsWith("LMS") || methodName.endsWith("RLS"))
                model.numOfFilterTrainingSymbols = 100;

              static if(methodName.split("_")[0].endsWith("FHF"))
                dirs ~= "TXP%s_inr%s_%s_B%s_G%s".format(model.pa.TX_POWER, model.INR, methodName, model.basisFuncsSelection.imageMargin, model.basisFuncsSelection.noiseMargin);
              else
                dirs ~= "TXP%s_inr%s_%s_orth%s".format(model.pa.TX_POWER, model.INR, methodName, model.orthogonalizer.numOfTrainingSymbols);
            }

            foreach(i; 0 .. models.length)
                taskList.append((Model m, string dir){
                    JSONValue[] resList;

                    // 最初の一回は普通にやる
                    resList ~= mainImpl!methodName(m, dir);

                    static if(methodName.startsWith("SFHF") || methodName.startsWith("S1FHF") || methodName.startsWith("S2FHF"))
                    {
                        // writeln(mainImpl!methodName(m, dir)["training_symbols_per_second"]);
                        enum K = 10;    // 試行回数
                        uint sumOfSuccFreq;
                        JSONValue[] selectingRatioList;
                        foreach(j; 0 .. K){
                            m.rndSeed += 100;   // seed値を100ずつ足していく
                            auto res = mainImpl!methodName(m, null);
                            resList ~= res;
                            auto cnt = res["filterSpec"]["selectingIsSuccess"].array.map!(a => a.type == JSON_TYPE.TRUE ? 1 : 0).sum();
                            sumOfSuccFreq += cnt;
                            selectingRatioList ~= JSONValue(cnt / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling));
                        }

                        JSONValue jv = cast(JSONValue[string])null;

                        jv["cancellations"] = resList.map!(a => a["cancellation_dB"]).array();
                        
                        jv["selecting"] = (){
                            JSONValue[string] vv;
                            vv["selectingRatio"] = sumOfSuccFreq / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * K);
                            vv["selectingRatioList"] = selectingRatioList;

                            size_t failedCNT, failedCNTPLQ;
                            size_t matchCNT, matchCNTPLQ;
                            size_t overCNT, overCNTPLQ;
                            real idealCostRLS = 0;
                            real idealCostLMS = 0;
                            real actualCostRLS = 0;
                            real actualCostLMS = 0;
                            size_t distDim = 0;
                            foreach(k, ref ev; resList){
                                auto reqs = ev["filterSpec"]["actualRequiredBasisFuncs"].array;
                                auto used = ev["filterSpec"]["selectedBasisFuncs"].array;
                                foreach(f; 0 .. m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling){
                                    auto r = reqs[f].array;
                                    auto u = used[f].array;

                                    if(distDim == 0) distDim = r.length;

                                    // idealCostに追加する
                                    size_t reqsP, usedP;
                                    foreach(p; 0 .. r.length){
                                        if(r[p].type == JSON_TYPE.TRUE) ++reqsP;
                                        if(u[p].type == JSON_TYPE.TRUE) ++usedP;
                                    }
                                    idealCostRLS += 4*(reqsP^^2) + 4*reqsP;
                                    idealCostLMS += 2*reqsP;
                                    actualCostRLS += 4*(usedP^^2) + 4*usedP;
                                    actualCostLMS += 2*usedP;

                                    // 推定ミスした周波数のカウント
                                    foreach(p; 0 .. r.length) if(r[p].type == JSON_TYPE.TRUE && u[p].type == JSON_TYPE.FALSE) { ++failedCNT; break; }

                                    // 過大評価した周波数のカウント
                                    foreach(p; 0 .. r.length) if(r[p].type == JSON_TYPE.FALSE && u[p].type == JSON_TYPE.TRUE){ ++overCNT; break; }
                                    
                                    // 完全推定できた周波数のカウント
                                    foreach(p; 0 .. r.length){
                                        if(r[p].type != u[p].type) break;
                                        if(p == r.length - 1) ++matchCNT;
                                    }

                                    // p>qについて
                                    // 推定ミスした周波数のカウント
                                    foreach(p; 0 .. r.length) if(CompleteDistorter!().indexOfConjugated(p) >= p) if(r[p].type == JSON_TYPE.TRUE && u[p].type == JSON_TYPE.FALSE) { ++failedCNTPLQ; break; }

                                    // 過大評価した周波数のカウント
                                    foreach(p; 0 .. r.length) if(CompleteDistorter!().indexOfConjugated(p) >= p) if(r[p].type == JSON_TYPE.FALSE && u[p].type == JSON_TYPE.TRUE){ ++overCNTPLQ; break; }
                                    
                                    // 完全推定できた周波数のカウント
                                    bool bBreaked = false;
                                    foreach(p; 0 .. r.length) if(CompleteDistorter!().indexOfConjugated(p) >= p) {
                                        if(r[p].type != u[p].type){
                                            bBreaked = true;
                                            break;
                                        }
                                    }
                                    if(!bBreaked) ++matchCNTPLQ;
                                }
                            }
                            vv["failedRatio"] = failedCNT / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);
                            vv["matchRatio"] = matchCNT / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);
                            vv["overRatio"] = overCNT / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);

                            vv["failedPLQRatio"] = failedCNTPLQ / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);
                            vv["matchPLQRatio"] = matchCNTPLQ / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);
                            vv["overPLQRatio"] = overCNTPLQ / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);

                            immutable size_t nfft = m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling;
                            immutable size_t nsym = (m.ofdm.numOfFFT + m.ofdm.numOfCP) * m.ofdm.scaleOfUpSampling;
                            vv["idealCostRLS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + (idealCostRLS / resList.length))/nsym;
                            vv["idealCostLMS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + (idealCostLMS / resList.length))/nsym;
                            vv["actualCostRLS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + (actualCostRLS / resList.length))/nsym;
                            vv["actualCostLMS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + (actualCostLMS / resList.length))/nsym;
                            vv["nonSelCostRLS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + 4*(distDim^^2 + distDim) * nfft)/nsym;
                            vv["nonSelCostLMS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + 2 * distDim * nfft)/nsym;

                            return vv;
                        }();
                        
                        auto file = File(buildPath(dir, "allResult.json"), "w");
                        file.write(jv.toPrettyString(JSONOptions.specialFloatLiterals));
                    }

                }, models[i], buildPath("results", dirs[i]));
            // foreach(i; 0 .. models.length)
            //     mainImpl!methodName(models[i], buildPath("results", dirs[i]));
        }

    //auto scheduler = new MPITaskScheduler();
    JobEnvironment env;
    pushArrayJob(taskList, env);
}


void main()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    // jobRun(1, 0, {
        mainJob();
    // });
}
