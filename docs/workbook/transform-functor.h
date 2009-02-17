// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


class TranslateTransform : public TransformBase<TranslateTransform> {
  double m_xtrans, m_ytrans;
public:
  TranslateTransform(double x_translation, double y_translation) : 
    m_xtrans( x_translation ) , m_ytrans( y_translation ) {}
  
  // Given a pixel coordinate in the ouput image, return
  // a pixel coordinate in the input image.
  inline Vector2 reverse(const Vector2 &p) const {
    return Vector2( p(0) - m_xtrans, p(1) - m_ytrans );
  }
  
  // Given a pixel coordinate in the input image, return
  // a pixel coordinate in the output image.
  inline Vector2 forward(const Vector2 &p) const {
    return Vector2( p(0) + m_xtrans, p(1) + m_ytrans );
  }
};
