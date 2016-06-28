#include <itpp/itcomm.h>
#include <cstdint>

using namespace itpp;

#ifdef _WIN64
   //define something for Windows (64-bit)
#elif _WIN32
   //define something for Windows (32-bit)
#elif __APPLE__
    #include "TargetConditionals.h"
    #if TARGET_OS_IPHONE && TARGET_IPHONE_SIMULATOR
        // define something for simulator   
    #elif TARGET_OS_IPHONE
        // define something for iphone  
    #else
        #define TARGET_OS_OSX 1
        #define THREAD_LOCAL __thread
    #endif
#elif __linux
    #define THREAD_LOCAL thread_local
#elif __unix // all unices not caught above
    // Unix
#elif __posix
    // POSIX
#endif


namespace WITPP
{
  itpp::bvec& to_bvec(itpp::bvec& dst, unsigned char const * src, uint32_t len)
  {
      dst.set_length(len);

      for(uint32_t i = 0; i < len; ++i)
          dst[i] = (src[i] != 0 ? 1 : 0);

      return dst;
  }


  template <typename F1, typename F2>
  itpp::Vec<std::complex<F1>>& to_cvec(itpp::Vec<std::complex<F1>>& dst, F2 const * src, uint32_t len)
  {
      dst.set_length(len);

      for(uint32_t i = 0; i < len; ++i){
          dst[i].real(src[i*2+0]);
          dst[i].imag(src[i*2+1]);
      }

      return dst;
  }



  void bvec_to_buf(itpp::bvec const & src, unsigned char * dst, uint32_t * len)
  {
      if(*len < src.length()){
          *len = 0;
          return;
      }

      *len = src.length();
      auto l = *len;
      for(uint32_t i = 0; i < l; ++i)
          dst[i] = src[i].value();
  }


  template <typename F>
  void cvec_to_buf(itpp::cvec const & src, F * dst, uint32_t * len)
  {
      if(*len < src.length()){
          *len = 0;
          return;
      }

      *len = src.length();
      auto l = *len;
      for(uint32_t i = 0; i < l; ++i){
          dst[i*2+0] = src[i].real();
          dst[i*2+1] = src[i].imag();
      }
  }

  /*
  namespace vec
  {
    void* makeObject() { return new itpp::vec(); }


    itpp::vec* makeObject(const double *ptr, const int size)
    {
        return new itpp::vec(ptr, size);
    }


    int length(const itpp::vec* v) { return v->length(); }
    void length(itpp::vec* v, int size) { v->set_length(size, true); }
    double* 
  }
  */


  namespace QAM
  {
    itpp::QAM* newObject(uint32_t mary) { return new itpp::QAM(mary); }
    void deleteObject(itpp::QAM* obj) { delete obj; }


    /**
    QAM変調された信号signalCpxを復調する．

    復調結果は一旦内部に保存され，その結果を得たい場合にはbitsとbits_lenに配列とその長さへのポインタを渡す．
    もし，復調結果を格納するために十分な容量がなければ，復調結果はbitsには格納されない．

    返却値：　現在内部で保持している復調後のビット列の長さ
    */
    uint32_t demodulate_bits(const itpp::QAM * obj, float const * signalCpx, uint32_t signal_len, unsigned char * bits, uint32_t * bits_len)
    {
        // auto& cv1 = tls<itpp::cvec, __LINE__, 0>;
        THREAD_LOCAL static itpp::cvec* _cv1 = nullptr;
        // auto& bv1 = tls<itpp::bvec, __LINE__, 0>;
        THREAD_LOCAL static itpp::bvec* _bv1 = nullptr;

        if(_cv1 == nullptr || _bv1 == nullptr)
        {
            _cv1 = new itpp::cvec();
            _bv1 = new itpp::bvec();
        }

        itpp::cvec& cv1 = *_cv1;
        itpp::bvec& bv1 = *_bv1;

        // 復調処理
        if(obj != nullptr && signalCpx != nullptr && signal_len != 0)
        {
            cv1 = to_cvec(cv1, signalCpx, signal_len);
            obj->demodulate_bits(cv1, bv1);
        }

        // 内部バッファをメモリに書き込む
        if((obj == nullptr || signalCpx == nullptr || signal_len == 0)
          && (bits != nullptr && bits_len != nullptr)
          && (*bits_len >= bv1.length()))
        {
            bvec_to_buf(bv1, bits, bits_len);
        }

        return bv1.length();
    }


    /**
    bit列をQAM変調する．

    変調結果は一旦内部に保存され，その結果を得たい場合にはsigsとsigs_lenに配列とその長さへのポインタを渡す．
    もし，変調結果を格納するために十分な容量がなければ，変調結果はbitsには格納されない．

    返却値：　現在内部で保持している変調後の信号の長さ
    */
    uint32_t modulate_bits(const itpp::QAM * obj, unsigned char const * bits, uint32_t bits_len, float * sigs, uint32_t * sigs_len)
    {
        // auto& cv1 = tls<itpp::cvec, __LINE__, 0>;
        THREAD_LOCAL static itpp::cvec* _cv1;
        // auto& bv1 = tls<itpp::bvec, __LINE__, 0>;
        THREAD_LOCAL static itpp::bvec* _bv1;

        if(_cv1 == nullptr || _bv1 == nullptr)
        {
            _cv1 = new itpp::cvec();
            _bv1 = new itpp::bvec();
        }

        itpp::cvec& cv1 = *_cv1;
        itpp::bvec& bv1 = *_bv1;

        // 変調処理
        if(obj != nullptr && bits != nullptr && bits_len != 0)
        {
            bv1 = to_bvec(bv1, bits, bits_len);
            obj->modulate_bits(bv1, cv1);
        }

        // 内部バッファをメモリに書き込む
        if((obj == nullptr || bits == nullptr || bits_len == 0)
          && (sigs != nullptr && sigs_len != nullptr)
          && (*sigs_len >= cv1.length()))
        {
            cvec_to_buf(cv1, sigs, sigs_len);
        }

        return cv1.length();
    }
  }


  namespace OFDM
  {
    itpp::OFDM* newObject(uint32_t inNfft, uint32_t inNcp, uint32_t inNupsample)
    {
        return new itpp::OFDM(inNfft, inNcp, inNupsample);
    }

    void deleteObject(itpp::OFDM* obj) { delete obj; }

    uint32_t modulate(itpp::OFDM * obj, float const * input, uint32_t inp_len, float * output, uint32_t * out_len)
    {
        THREAD_LOCAL static itpp::cvec* _cv1;
        THREAD_LOCAL static itpp::cvec* _cv2;

        if(_cv1 == nullptr || _cv2 == nullptr)
        {
            _cv1 = new itpp::cvec();
            _cv2 = new itpp::cvec();
        }

        itpp::cvec& cv1 = *_cv1;
        itpp::cvec& cv2 = *_cv2;

        if(obj != nullptr && input != nullptr && inp_len != 0)
        {
            cv1 = to_cvec(cv1, input, inp_len);
            obj->modulate(cv1, cv2);
        }

        if((obj == nullptr || input == nullptr || inp_len == 0)
          && (output != nullptr && out_len != nullptr)
          && (*out_len >= cv2.length()))
        {
            cvec_to_buf(cv2, output, out_len);
        }

        return cv2.length();
    }


    uint32_t demodulate(itpp::OFDM * obj, float const * input, uint32_t inp_len, float * output, uint32_t * out_len)
    {
        THREAD_LOCAL static itpp::cvec* _cv1;
        THREAD_LOCAL static itpp::cvec* _cv2;

        if(_cv1 == nullptr || _cv2 == nullptr)
        {
            _cv1 = new itpp::cvec();
            _cv2 = new itpp::cvec();
        }

        itpp::cvec& cv1 = *_cv1;
        itpp::cvec& cv2 = *_cv2;

        if(obj != nullptr && input != nullptr && inp_len != 0)
        {
            cv1 = to_cvec(cv1, input, inp_len);
            obj->demodulate(cv1, cv2);
        }

        if((obj == nullptr || input == nullptr || inp_len == 0)
          && (output != nullptr && out_len != nullptr)
          && (*out_len >= cv2.length()))
        {
            cvec_to_buf(cv2, output, out_len);
        }

        return cv2.length();
    }
  }
}