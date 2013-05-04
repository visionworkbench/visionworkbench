// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
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
