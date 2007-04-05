#include <vw/Core/ProgressCallback.h>

namespace vw {
  namespace {
    ProgressCallback g_dummy_progress_callback_instance;
  }

  const ProgressCallback &ProgressCallback::dummy_instance() { return g_dummy_progress_callback_instance; }
}
