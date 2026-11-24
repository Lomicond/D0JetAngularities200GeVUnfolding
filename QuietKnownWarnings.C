#include "TError.h"
#include <cstring>

static ErrorHandlerFunc_t previousKnownWarningHandler = 0;

void KnownWarningHandler(
    Int_t level, Bool_t abortFlag,
    const char *location, const char *message)
{
    if (level == kWarning && !abortFlag && location && message) {
        const Bool_t knownLocation =
            std::strcmp(location, "TROOT::Append") == 0;

        const Bool_t knownMessage =
            std::strcmp(message,
                "Replacing existing TH1: h_total (Potential memory leak).") == 0 ||
            std::strcmp(message,
                "Replacing existing TH1: h_passed (Potential memory leak).") == 0;

        if (knownLocation && knownMessage)
            return;
    }

    if (previousKnownWarningHandler)
        previousKnownWarningHandler(level, abortFlag, location, message);
    else
        DefaultErrorHandler(level, abortFlag, location, message);
}

void QuietKnownWarnings()
{
    if (GetErrorHandler() != KnownWarningHandler)
        previousKnownWarningHandler = SetErrorHandler(KnownWarningHandler);
}
