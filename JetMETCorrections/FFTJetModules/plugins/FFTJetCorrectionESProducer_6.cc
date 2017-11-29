#include "JetMETCorrections/FFTJetModules/plugins/FFTJetCorrectionESProducer.h"

//
// define this as a plug-in
//
typedef FFTJetCorrectionESProducer<fftcorrtypes::CHS0Sys> FFTCHS0SysCorrectionESProducer;
typedef FFTJetCorrectionESProducer<fftcorrtypes::CHS1Sys> FFTCHS1SysCorrectionESProducer;
typedef FFTJetCorrectionESProducer<fftcorrtypes::CHS2Sys> FFTCHS2SysCorrectionESProducer;
typedef FFTJetCorrectionESProducer<fftcorrtypes::CHS3Sys> FFTCHS3SysCorrectionESProducer;
typedef FFTJetCorrectionESProducer<fftcorrtypes::CHS4Sys> FFTCHS4SysCorrectionESProducer;
typedef FFTJetCorrectionESProducer<fftcorrtypes::CHS5Sys> FFTCHS5SysCorrectionESProducer;
typedef FFTJetCorrectionESProducer<fftcorrtypes::CHS6Sys> FFTCHS6SysCorrectionESProducer;
typedef FFTJetCorrectionESProducer<fftcorrtypes::CHS7Sys> FFTCHS7SysCorrectionESProducer;
typedef FFTJetCorrectionESProducer<fftcorrtypes::CHS8Sys> FFTCHS8SysCorrectionESProducer;
typedef FFTJetCorrectionESProducer<fftcorrtypes::CHS9Sys> FFTCHS9SysCorrectionESProducer;

DEFINE_FWK_EVENTSETUP_MODULE(FFTCHS0SysCorrectionESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(FFTCHS1SysCorrectionESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(FFTCHS2SysCorrectionESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(FFTCHS3SysCorrectionESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(FFTCHS4SysCorrectionESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(FFTCHS5SysCorrectionESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(FFTCHS6SysCorrectionESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(FFTCHS7SysCorrectionESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(FFTCHS8SysCorrectionESProducer);
DEFINE_FWK_EVENTSETUP_MODULE(FFTCHS9SysCorrectionESProducer);
