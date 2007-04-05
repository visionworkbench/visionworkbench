#ifndef PROGRESSCALLBACK_H
#define PROGRESSCALLBACK_H

// System includes
#include <math.h>

// Design considerations:
// Should halt throw an exception instead?

namespace vw {
  class ProgressCallback {
  protected:
    bool m_abort_requested;
  public:
    ProgressCallback() : m_abort_requested( false ) {}

    // Reporting functions
    // Subclasses should reimplement where appropriate
    //
    // progress is from 0 (not done) to 1 (finished)
    // 
    virtual void report_progress(double /*progress*/) const {}
    virtual void report_aborted() const {}
    virtual void report_finished() const {}

    // Helper method which computes progress and calls report_progress
    void report_fractional_progress(double n, double total) const {
      report_progress(fabs(total) < 1e-30 ? 0 : n/total);
    }

    // Has an abort been requested?
    virtual bool abort_requested() const { return m_abort_requested; }

    // Request abort
    virtual void request_abort() { m_abort_requested = true; }
    virtual ~ProgressCallback() {}
    static const ProgressCallback &dummy_instance();
  };

  class SubProgressCallback : public ProgressCallback {
  protected:
    const ProgressCallback &m_parent;
    const double m_from;
    const double m_to;
  public:
    SubProgressCallback(const ProgressCallback &parent,
                        double from, double to) :
      m_parent(parent), m_from(from), m_to(to) {}
    virtual void report_progress(double progress) const {
      double parent_progress = m_from + (m_to - m_from)*progress;
      m_parent.report_progress(parent_progress);
    }
    virtual bool abort_requested() const { return m_parent.abort_requested(); }
    virtual ~SubProgressCallback() {}
  };
}

#endif
