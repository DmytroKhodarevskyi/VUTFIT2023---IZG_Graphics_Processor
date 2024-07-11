/*!
 * @file
 * @brief This file contains implementation of gpu
 *
 * @author Tomáš Milet, imilet@fit.vutbr.cz
 */

#include <student/gpu.hpp>

#include <iostream>
#include <string>

struct Triangle{
  OutVertex points[3];
};

struct point{
  float x;
  float y;
  float z;
};

struct clipP{
  float x;
  float y;
  float z;
  float w;
  bool isOutside;
  bool isInside;
  Attribute attributes[maxAttributes];
};

struct ClippedTriangle{
  Triangle tr1;
  Triangle tr2;
  bool isClipped;
};

//MACRO
void p(std::string s) {std::cerr << s << std::endl;}
//MACRO

bool pixel_in_triangle_bar(float x, float y, Triangle triangle, float &lambdaA, float &lambdaB, float &lambdaC){
  // vypočítat barycentrické souřadnice
  // pokud jsou všechny >= 0, je bod uvnitř trojúhelníku
  // pokud je nějaká < 0, bod je vně trojúhelníku
  // pokud je nějaká == 0, bod leží na hraně trojúhelníku

  point pixel = {x + 0.5f, y + 0.5f}; // střed pixelu

  point p1 = {triangle.points[0].gl_Position.x, triangle.points[0].gl_Position.y, triangle.points[0].gl_Position.z};
  point p2 = {triangle.points[1].gl_Position.x, triangle.points[1].gl_Position.y, triangle.points[1].gl_Position.z};
  point p3 = {triangle.points[2].gl_Position.x, triangle.points[2].gl_Position.y, triangle.points[2].gl_Position.z};

  glm::mat3 matrix = {p1.x, p2.x, p3.x,
                      p1.y, p2.y, p3.y,
                       1,    1,    1};

  glm::mat3 la = {pixel.x, p2.x, p3.x,
                  pixel.y, p2.y, p3.y,
                  1,       1,    1};

  glm::mat3 lb = {p1.x, pixel.x, p3.x,
                  p1.y, pixel.y, p3.y,
                  1,       1,    1};

  glm::mat3 lc = {p1.x, p2.x, pixel.x,
                  p1.y, p2.y, pixel.y,
                  1,       1,    1};


  float det = glm::determinant(matrix);
  lambdaA = glm::determinant(la);
  lambdaB = glm::determinant(lb);
  lambdaC = glm::determinant(lc);

  lambdaA = lambdaA/det;
  lambdaB = lambdaB/det;
  lambdaC = lambdaC/det;

  if (lambdaA >= 0. && lambdaB >= 0. && lambdaC >= 0.) {
    return true;
  }

  return false;
}


//UNUSED
bool pixel_in_triangle_edge(float x, float y, Triangle triangle, float &lambdaA, float &lambdaB, float &lambdaC){

  point A = {triangle.points[0].gl_Position.x, triangle.points[0].gl_Position.y, triangle.points[0].gl_Position.z};
  point B = {triangle.points[1].gl_Position.x, triangle.points[1].gl_Position.y, triangle.points[1].gl_Position.z};
  point C = {triangle.points[2].gl_Position.x, triangle.points[2].gl_Position.y, triangle.points[2].gl_Position.z};

  point P = {x + 0.5f, y + 0.5f, 0};

  point v0 = {C.x - A.x, C.y - A.y};
  point v1 = {B.x - A.x, B.y - A.y};
  point v2 = {P.x - A.x, P.y - A.y};

  // Calculate dot products
  float dot00 = v0.x * v0.x + v0.y * v0.y;
  float dot01 = v0.x * v1.x + v0.y * v1.y;
  float dot02 = v0.x * v2.x + v0.y * v2.y;
  float dot11 = v1.x * v1.x + v1.y * v1.y;
  float dot12 = v1.x * v2.x + v1.y * v2.y;

  // Calculate barycentric coordinates
  float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
  lambdaA = (dot11 * dot02 - dot01 * dot12) * invDenom;
  lambdaB = (dot00 * dot12 - dot01 * dot02) * invDenom;
  lambdaC = 1.0f - lambdaA - lambdaB;

  // Check if the point is inside the triangle
  return lambdaA >= 0.0f && lambdaB >= 0.0f && lambdaC >= 0.0f;
}
//UNUSED


bool backface_check(Triangle tr){

  point A = {tr.points[0].gl_Position.x, tr.points[0].gl_Position.y, tr.points[0].gl_Position.z};
  point B = {tr.points[1].gl_Position.x, tr.points[1].gl_Position.y, tr.points[1].gl_Position.z};
  point C = {tr.points[2].gl_Position.x, tr.points[2].gl_Position.y, tr.points[2].gl_Position.z};

  glm::vec3 AB = {B.x - A.x, B.y - A.y, B.z - A.z};
  glm::vec3 AC = {C.x - A.x, C.y - A.y, C.z - A.z};

  auto cross = glm::cross(AB, AC);

  if(cross.z > 0){ //cw
    return true;
  }

  return false; //ccw
}

float depth_interpolation(float lambdaA, float lambdaB, float lambdaC, Triangle triangle){
  float z1 = triangle.points[0].gl_Position.z;
  float z2 = triangle.points[1].gl_Position.z;
  float z3 = triangle.points[2].gl_Position.z;

  float depth = lambdaA * z1 + lambdaB * z2 + lambdaC * z3;

  return depth;
}

void color_clamp(OutFragment &oF){
  oF.gl_FragColor.x = glm::max(0.0f, glm::min(1.0f, oF.gl_FragColor.x));
  oF.gl_FragColor.y = glm::max(0.0f, glm::min(1.0f, oF.gl_FragColor.y));
  oF.gl_FragColor.z = glm::max(0.0f, glm::min(1.0f, oF.gl_FragColor.z));
  oF.gl_FragColor.w = glm::max(0.0f, glm::min(1.0f, oF.gl_FragColor.w));
}

float fclamp(float color, float min, float max){
  return glm::max(min, glm::min(max, color));
}

void interpolate_colors(Triangle const& tr, AttributeType const* vs2fs, InFragment& inF, float lambdaA, float lambdaB, float lambdaC) {
    for (int i = 0; i < maxAttributes; i++) {
        AttributeType attributeType = vs2fs[i];

        // Get the attribute values of the triangle vertices
        auto trAattribute = tr.points[0].attributes[i];
        auto trBattribute = tr.points[1].attributes[i];
        auto trCattribute = tr.points[2].attributes[i];

        if (attributeType != AttributeType::EMPTY && attributeType != AttributeType::UINT &&
            attributeType != AttributeType::UVEC2 && attributeType != AttributeType::UVEC3 &&
            attributeType != AttributeType::UVEC4) {
            // Interpolate only non-empty and non-integer attributes

            float h1 = tr.points[0].gl_Position.w;
            float h2 = tr.points[1].gl_Position.w;
            float h3 = tr.points[2].gl_Position.w;

            float l1 = lambdaA;
            float l2 = lambdaB;
            float l3 = lambdaC;

            float s = l1/h1 + l2/h2 + l3/h3;

            float l13d = l1/(h1*s);
            float l23d = l2/(h2*s);
            float l33d = l3/(h3*s);

            if (attributeType == AttributeType::VEC2) {
                inF.attributes[i].v2 = l13d * trAattribute.v2 + l23d * trBattribute.v2 + l33d * trCattribute.v2;
            } else if (attributeType == AttributeType::VEC3) {
                inF.attributes[i].v3 = l13d * trAattribute.v3 + l23d * trBattribute.v3 + l33d * trCattribute.v3;
            } else if (attributeType == AttributeType::VEC4) {
                inF.attributes[i].v4 = l13d * trAattribute.v4 + l23d * trBattribute.v4 + l33d * trCattribute.v4;
            }  else if (attributeType == AttributeType::FLOAT) {
                inF.attributes[i].v1 = l13d * trAattribute.v1 + l23d * trBattribute.v1 + l33d * trCattribute.v1;
            }

            // inF.attributes[i].v1 = l13d * trAattribute.v1 + l23d * trBattribute.v1 + l33d * trCattribute.v1;
            // inF.attributes[i].v2 = l13d * trAattribute.v2 + l23d * trBattribute.v2 + l33d * trCattribute.v2;
            // inF.attributes[i].v3 = l13d * trAattribute.v3 + l23d * trBattribute.v3 + l33d * trCattribute.v3;
            // inF.attributes[i].v4 = l13d * trAattribute.v4 + l23d * trBattribute.v4 + l33d * trCattribute.v4;
        } else {
            // For integer attributes, use zero vertex values (provoking vertex)
            if (attributeType == AttributeType::EMPTY)
              continue;
            if (attributeType == AttributeType::UINT) {
                inF.attributes[i].v1 = trAattribute.v1;
                // inF.attributes[i].u1 = trAattribute.v1;
                inF.attributes[i].u1 = trAattribute.u1;
            } else if (attributeType == AttributeType::UVEC2) {
                inF.attributes[i].v2 = trAattribute.v2;
                // inF.attributes[i].u2 = trAattribute.v2;
                inF.attributes[i].u2 = trAattribute.u2;

            } else if (attributeType == AttributeType::UVEC3) {
                inF.attributes[i].v3 = trAattribute.v3;
                // inF.attributes[i].u3 = trAattribute.v3;
                inF.attributes[i].u3 = trAattribute.u3;

            } else if (attributeType == AttributeType::UVEC4) {
                inF.attributes[i].v4 = trAattribute.v4;
                // inF.attributes[i].u4 = trAattribute.v4;
                inF.attributes[i].u4 = trAattribute.u4;

            }
        }
    }
}

void PFO(Frame &frame, OutFragment &oF, float depth, uint32_t x, uint32_t y){
  uint32_t w = frame.width;
  uint32_t h = frame.height;

        float xcol = oF.gl_FragColor.x;
        float ycol = oF.gl_FragColor.y;
        float zcol = oF.gl_FragColor.z;
        float alpha = oF.gl_FragColor.a;

        if (depth < frame.depth[y * w + x]){

          if (alpha > 0.5f)
            frame.depth[y * w + x] = depth;

          uint8_t* bufColX = &(frame.color[(y * w + x) * 4]);
          uint8_t* bufColY = &(frame.color[(y * w + x) * 4 + 1]);
          uint8_t* bufColZ = &(frame.color[(y * w + x) * 4 + 2]);
          // uint8_t* bufColA = &(frame.color[(y * w + x) * 4 + 3]);

          *bufColX = fclamp((*bufColX/255.f)*(1-alpha) + (xcol)*(alpha),0.f,1.f)*255.f;
          *bufColY = fclamp((*bufColY/255.f)*(1-alpha) + (ycol)*(alpha),0.f,1.f)*255.f;
          *bufColZ = std::lroundf(fclamp((*bufColZ/255.f)*(1-alpha) + (zcol)*(alpha),0.f,1.f)*255.f);
          // *bufColA = alpha*255.f;

        }
}

void rasterize(Frame&frame, Triangle const&triangle, Program const&program,
               GPUMemory &mem, bool backFaceCulling){
      // spočítat hranice, trojúhelníka
      //
      uint32_t w = frame.width;
      uint32_t h = frame.height;
      
      if(backFaceCulling){
        if (!backface_check(triangle)){
          return;
        }
      }

      float lA, lB, lC;

      glm::vec4 A = triangle.points[0].gl_Position;
      glm::vec4 B = triangle.points[1].gl_Position;
      glm::vec4 C = triangle.points[2].gl_Position;

      float minX = std::floor(glm::min(glm::min(A.x, B.x), C.x));
      float maxX = std::ceil(glm::max(glm::max(A.x, B.x), C.x));
      float minY = std::floor(glm::min(glm::min(A.y, B.y), C.y));
      float maxY = std::ceil(glm::max(glm::max(A.y, B.y), C.y));

      if (minX < 0)
        minX = 0;

      if (minY < 0)
        minY = 0;

      if (maxX > w)
        maxX = w;

      if (maxY > h)
        maxY = h;

      // std::cerr << "begin rasterize..." << std::endl;

      for(float y = minY; y < maxY; y++){ // smyčka přes pixely
        for(float x = minX; x < maxX; x++){
          if(pixel_in_triangle_bar(x, y, triangle, lA, lB, lC)){

            InFragment inFragment;

            

            inFragment.gl_FragCoord = {x + 0.5f, y + 0.5f, 0, 1};
            inFragment.gl_FragCoord.z = depth_interpolation(lA, lB, lC, triangle);
            interpolate_colors(triangle, program.vs2fs, inFragment, lA, lB, lC);

            // std::cerr << "aaaaa " << inFragment.attributes[3].u1 << std::endl;

            // std::cerr << "normal " << inFragment.attributes[1].v3.x << " " 
            // << inFragment.attributes[1].v3.y << " " 
            // << inFragment.attributes[1].v3.z << std::endl;

            //  std::cerr << "position " << inFragment.attributes[0].v3.x << " " 
            //  << inFragment.attributes[0].v3.y << " " 
            //  << inFragment.attributes[0].v3.z << std::endl;

            OutFragment outFragment;
            ShaderInterface si;

            si.uniforms = mem.uniforms;
            si.textures = mem.textures;
            // p("pls1");


            program.fragmentShader(outFragment, inFragment, si);

            // p("pls");
            color_clamp(outFragment);
            PFO(frame, outFragment, inFragment.gl_FragCoord.z, x, y);
          }
        }
      }
      // p("rasterize end");
}

void clear(GPUMemory &mem,ClearCommand cmd){

  if (cmd.clearColor){
    for (uint32_t i = 0; i < mem.framebuffer.width*mem.framebuffer.height; ++i){
      mem.framebuffer.color[i*4+0] = cmd.color.r * 255.f;
      mem.framebuffer.color[i*4+1] = cmd.color.g * 255.f;
      mem.framebuffer.color[i*4+2] = cmd.color.b * 255.f;
      mem.framebuffer.color[i*4+3] = cmd.color.a * 255.f;
    }
  }

  if (cmd.clearDepth){
    for (uint32_t i = 0; i < mem.framebuffer.width*mem.framebuffer.height; ++i){
      mem.framebuffer.depth[i] = cmd.depth;
    }
  }
  
}

void readAttributes(GPUMemory const&mem, DrawCommand cmd, uint32_t j, InVertex &inVertex)
{
  for (int k = 0; k < maxAttributes; ++k){

            if (cmd.vao.vertexAttrib[k].bufferID < 0){
              // std::cout << "buffer is nullptr" << std::endl;
              continue;
            }

        Buffer attribsBuffer = mem.buffers[cmd.vao.vertexAttrib[k].bufferID];
        int stride = cmd.vao.vertexAttrib[k].stride;
        int offset = cmd.vao.vertexAttrib[k].offset;
        AttributeType type = cmd.vao.vertexAttrib[k].type;

        uint8_t* dataPtr = reinterpret_cast<uint8_t*>(const_cast<void*>(attribsBuffer.data));

          if (type == AttributeType::FLOAT){
            // p("float");
            inVertex.attributes[k].v1 = (float) ((float*) (dataPtr + offset + j*stride))[0];
          }

          if (type == AttributeType::VEC2){
            // p("vec2");
            inVertex.attributes[k].v2 = (glm::vec2) ((glm::vec2*) (dataPtr + offset + j*stride))[0];
          }

          if (type == AttributeType::VEC3){
            // p("vec3");
            inVertex.attributes[k].v3 = (glm::vec3) ((glm::vec3*) (dataPtr + offset + j*stride))[0];
          }                                     

          if (type == AttributeType::VEC4){
            // p("vec4");
            inVertex.attributes[k].v4 = (glm::vec4) ((glm::vec4*) (dataPtr + offset + j*stride))[0];
          }

          if (type == AttributeType::UINT){
            // p("uint");
            inVertex.attributes[k].u1 = (uint32_t) ((uint32_t*) (dataPtr + offset + j*stride))[0];
          }

          if (type == AttributeType::UVEC2){
            // p("uvec2");
            inVertex.attributes[k].u2 = (glm::uvec2) ((glm::uvec2*) (dataPtr + offset + j*stride))[0];
          }

          if (type == AttributeType::UVEC3){
            // p("uvec3");
            inVertex.attributes[k].u3 = (glm::uvec3) ((glm::uvec3*) (dataPtr + offset + j*stride))[0];
          }

          if (type == AttributeType::UVEC4){
            // p("uvec4");
            inVertex.attributes[k].u4 = (glm::uvec4) ((glm::uvec4*) (dataPtr + offset + j*stride))[0];
          }
        }
}

void loadVertex(InVertex &inVertex, GPUMemory const&mem, int vertnum, DrawCommand cmd, uint32_t drawID){

  // DrawCommand cmd = cb.commands[i].data.drawCommand;

  if (cmd.vao.indexBufferID < 0 )
    {
        inVertex.gl_VertexID = vertnum;

        readAttributes(mem, cmd, inVertex.gl_VertexID, inVertex);

        inVertex.gl_DrawID = drawID;
    }
    else
    {
        OutVertex outVertex;

        if (cmd.vao.indexType == IndexType::UINT32){

          Buffer indexBuffer = mem.buffers[cmd.vao.indexBufferID];

          uint8_t* dataPtr = reinterpret_cast<uint8_t*>(const_cast<void*>(indexBuffer.data));


          uint32_t* index = (uint32_t*) (dataPtr);

          inVertex.gl_VertexID = (uint32_t) index[vertnum + cmd.vao.indexOffset/4];
        }

        if (cmd.vao.indexType == IndexType::UINT16){

          Buffer indexBuffer = mem.buffers[cmd.vao.indexBufferID];

          uint8_t* dataPtr = reinterpret_cast<uint8_t*>(const_cast<void*>(indexBuffer.data));

          uint16_t* index = (uint16_t*) (dataPtr);

          inVertex.gl_VertexID = (uint16_t) index[vertnum + cmd.vao.indexOffset/2];


        }

        if (cmd.vao.indexType == IndexType::UINT8){

          Buffer indexBuffer = mem.buffers[cmd.vao.indexBufferID];

          uint8_t* dataPtr = reinterpret_cast<uint8_t*>(const_cast<void*>(indexBuffer.data));

          uint8_t index = dataPtr[vertnum + cmd.vao.indexOffset];
          inVertex.gl_VertexID = index;
        }

        readAttributes(mem, cmd, inVertex.gl_VertexID, inVertex);

    }
}

void loadTriangle(Triangle&triangle, uint32_t tId, GPUMemory const&mem, DrawCommand &cmd, uint32_t drawID){
      for(int i = 0; i < 3; i++){ // smyčka přes vrcholy trojúhelníku
        InVertex inVertex;
        loadVertex(inVertex, mem, tId+i, cmd, drawID);
        VertexShader vs = mem.programs[cmd.programID].vertexShader;
        // p("im dead 2?");
        ShaderInterface si;

        si.uniforms = mem.uniforms;
        // si.uniforms[1] = mem.uniforms[1];
        vs(triangle.points[i], inVertex, si);
        // p("im dead 3?");
      }
    }

//UNUSED
void draw(GPUMemory &mem, int i, CommandBuffer &cb, int drawID){

    DrawCommand cmd = cb.commands[i].data.drawCommand;
    VertexShader vs = mem.programs[cmd.programID].vertexShader;

    if (cmd.vao.indexBufferID < 0 )
    {
      for(uint32_t j = 0; j < cmd.nofVertices; ++j){

        InVertex inVertex;
        OutVertex outVertex;
        inVertex.gl_VertexID = j;

        readAttributes(mem, cmd, inVertex.gl_VertexID, inVertex);

        inVertex.gl_DrawID = drawID;

        ShaderInterface si;
        vs(outVertex,inVertex,si);
      }
    }
    else
    {
      for(uint32_t j = 0; j < cmd.nofVertices; ++j){
        InVertex inVertex;
        OutVertex outVertex;

        if (cmd.vao.indexType == IndexType::UINT32){

          Buffer indexBuffer = mem.buffers[cmd.vao.indexBufferID];

          uint8_t* dataPtr = reinterpret_cast<uint8_t*>(const_cast<void*>(indexBuffer.data));

          uint32_t* index = (uint32_t*) (dataPtr);

          inVertex.gl_VertexID = (uint32_t) index[j + cmd.vao.indexOffset/4];
        }

        if (cmd.vao.indexType == IndexType::UINT16){

          Buffer indexBuffer = mem.buffers[cmd.vao.indexBufferID];

          uint8_t* dataPtr = reinterpret_cast<uint8_t*>(const_cast<void*>(indexBuffer.data));

          uint16_t* index = (uint16_t*) (dataPtr);

          inVertex.gl_VertexID = (uint16_t) index[j + cmd.vao.indexOffset/2];


        }

        if (cmd.vao.indexType == IndexType::UINT8){

          Buffer indexBuffer = mem.buffers[cmd.vao.indexBufferID];

          uint8_t* dataPtr = reinterpret_cast<uint8_t*>(const_cast<void*>(indexBuffer.data));

          uint8_t index = dataPtr[j];
          inVertex.gl_VertexID = index;
        }

        readAttributes(mem, cmd, inVertex.gl_VertexID, inVertex);

        ShaderInterface si;
        vs(outVertex,inVertex,si);
      }
    }
}
//UNUSED

void perspectiveDivision(Triangle&triangle) {

  for (int i = 0; i < 3; i++) {

    float Xpos = triangle.points[i].gl_Position.x;
    float Ypos = triangle.points[i].gl_Position.y;
    float Zpos = triangle.points[i].gl_Position.z;
    float Wpos = triangle.points[i].gl_Position.w;
    
    triangle.points[i].gl_Position.x = Xpos / Wpos;
    triangle.points[i].gl_Position.y = Ypos / Wpos;
    triangle.points[i].gl_Position.z = Zpos / Wpos;
  }
}

void viewportTrans(Triangle &triangle, float width, float height){

  for (int i = 0; i < 3; i++){
    triangle.points[i].gl_Position.x = (triangle.points[i].gl_Position.x + 1) * width / 2;
    triangle.points[i].gl_Position.y = (triangle.points[i].gl_Position.y + 1) * height / 2;
    // triangle.points[i].gl_Position.z = (triangle.points[i].gl_Position.z + 1) / 2;
    // triangle.points[i].gl_Position.w = (triangle.points[i].gl_Position.w + 1) / 2;
  }
}

void set_coords(Triangle &tr, int point, float x, float y, float z, float w){
  tr.points[point].gl_Position.x = x;
  tr.points[point].gl_Position.y = y;
  tr.points[point].gl_Position.z = z;
  tr.points[point].gl_Position.w = w;
}

void set_attribs(clipP A, ClippedTriangle &cltr, int num_tr, int point_num){
  if (num_tr == 1){
    for (int i = 0; i < maxAttributes; i++){
      cltr.tr1.points[point_num].attributes[i] = A.attributes[i];
    }
    // cltr.tr1.points[point_num].gl_Position = A.attributes[i];
  }
  else if (num_tr == 2){
    for (int i = 0; i < maxAttributes; i++){
      cltr.tr2.points[point_num].attributes[i] = A.attributes[i];
    }
    // cltr.tr2.points[point_num].gl_Position = A.attributes[i];
  }
  // for (int i = 0; i < maxAttributes; i++){
    // cltr.tr1.points[point_num].attributes[i] = A.attributes[i];
  // }
}

void interpolate_attribs (clipP inside, clipP outside, ClippedTriangle &cltr,
                        int num_tr, int point_num, float t) {
  for (int i = 0; i < maxAttributes; i++){
  //         cltr.tr1.points[0].attributes[i] = A.attributes[i];
  if (num_tr == 1){
      float ins_atr_v1 = inside.attributes[i].v1;
      float out_atr_v1 = outside.attributes[i].v1;
      cltr.tr1.points[point_num].attributes[i].v1 = ins_atr_v1 + t * (out_atr_v1 - ins_atr_v1);
      glm::vec2 ins_atr_v2 = inside.attributes[i].v2;
      glm::vec2 out_atr_v2 = outside.attributes[i].v2;
      cltr.tr1.points[point_num].attributes[i].v2 = ins_atr_v2 + t * (out_atr_v2 - ins_atr_v2);
      glm::vec3 ins_atr_v3 = inside.attributes[i].v3;
      glm::vec3 out_atr_v3 = outside.attributes[i].v3;
      cltr.tr1.points[point_num].attributes[i].v3 = ins_atr_v3 + t * (out_atr_v3 - ins_atr_v3);
      glm::vec4 ins_atr_v4 = inside.attributes[i].v4;
      glm::vec4 out_atr_v4 = outside.attributes[i].v4;
      cltr.tr1.points[point_num].attributes[i].v4 = ins_atr_v4 + t * (out_atr_v4 - ins_atr_v4);

      cltr.tr1.points[point_num].attributes[i].u1 = inside.attributes[i].u1;
  }
  if (num_tr == 2){
      float ins_atr_v1 = inside.attributes[i].v1;
      float out_atr_v1 = outside.attributes[i].v1;
      cltr.tr2.points[point_num].attributes[i].v1 = ins_atr_v1 + t * (out_atr_v1 - ins_atr_v1);
      glm::vec2 ins_atr_v2 = inside.attributes[i].v2;
      glm::vec2 out_atr_v2 = outside.attributes[i].v2;
      cltr.tr2.points[point_num].attributes[i].v2 = ins_atr_v2 + t * (out_atr_v2 - ins_atr_v2);
      glm::vec3 ins_atr_v3 = inside.attributes[i].v3;
      glm::vec3 out_atr_v3 = outside.attributes[i].v3;
      cltr.tr2.points[point_num].attributes[i].v3 = ins_atr_v3 + t * (out_atr_v3 - ins_atr_v3);
      glm::vec4 ins_atr_v4 = inside.attributes[i].v4;
      glm::vec4 out_atr_v4 = outside.attributes[i].v4;
      cltr.tr2.points[point_num].attributes[i].v4 = ins_atr_v4 + t * (out_atr_v4 - ins_atr_v4);

      cltr.tr2.points[point_num].attributes[i].u1 = inside.attributes[i].u1;
  }
  // if (cltr,)
  }
}

void clipping(GPUMemory &mem, ClippedTriangle &cltr, Triangle &tr, int plane, Frame &frame){

    int inside = 0;
    int outside = 0;
    int onPlane = 0;

    clipP A, B, C;

    A.x = tr.points[0].gl_Position.x;
    A.y = tr.points[0].gl_Position.y;
    A.z = tr.points[0].gl_Position.z;
    A.w = tr.points[0].gl_Position.w;
    // A.attributes = tr.points[0].attributes.;
    for (int i = 0; i < maxAttributes; i++){
      A.attributes[i] = tr.points[0].attributes[i];
      // std::cerr << A.attributes[i].u1 << std::endl;   
    }

    // std::cerr << A.attributes[3].u1 << std::endl;


    B.x = tr.points[1].gl_Position.x;
    B.y = tr.points[1].gl_Position.y;
    B.z = tr.points[1].gl_Position.z;
    B.w = tr.points[1].gl_Position.w;
    for (int i = 0; i < maxAttributes; i++){
      B.attributes[i] = tr.points[1].attributes[i];
    }

    C.x = tr.points[2].gl_Position.x;
    C.y = tr.points[2].gl_Position.y;
    C.z = tr.points[2].gl_Position.z;
    C.w = tr.points[2].gl_Position.w;
    for (int i = 0; i < maxAttributes; i++){
      C.attributes[i] = tr.points[2].attributes[i];
    }

    // −Aw≤Az

    // glm::vec4 ClipPoint1, ClipPoint2;
    clipP ClipPoint1, ClipPoint2;


    if (-A.w <= A.z){
      A.isOutside = true;
      A.isInside = false;
      outside++;
    } else {
      A.isOutside = false;
      A.isInside = true;
      inside++;
    }

    if (-B.w <= B.z){
      B.isOutside = true;
      B.isInside = false;
      outside++;
    } else {
      B.isOutside = false;
      B.isInside = true;
      inside++;
    }

    if (-C.w <= C.z){
      C.isOutside = true;
      C.isInside = false;
      outside++;
    } else {
      C.isOutside = false;
      C.isInside = true;
      inside++;
    }

    if (inside == 3){
      return;
    }

    if (outside == 3){
      cltr.tr1 = tr;
    }

    if (outside == 1 && inside == 2){
      // if B is outside
      // (-Aw - Az) / (Bw - Aw + Bz - Az) = t
      // ClipPoint1 = A + t(B - A)
      // (-Cw - Cz) / (Cw - Cw + Bz - Cz) = t2
      // ClipPoint2 = C + t2(B - C)
      
      // if C is outside
      // (-Aw - Az) / (Cw - Aw + Cz - Az) = t
      // ClipPoint1 = A + t(C - A)
      // (-Bw - Bz) / (Bw - Bw + Cz - Bz) = t2
      // ClipPoint2 = B + t2(C - B)

      // if A is outside
      // (-Bw - Bz) / (Aw - Bw + Az - Bz) = t
      // ClipPoint1 = B + t(A - B)
      // (-Cw - Cz) / (Cw - Cw + Az - Cz) = t2
      // ClipPoint2 = C + t2(A - C)


      if (A.isOutside)
      {
        float t = (-B.w - B.z) / (A.w - B.w + A.z - B.z);
        ClipPoint1.x = B.x + t * (A.x - B.x);
        ClipPoint1.y = B.y + t * (A.y - B.y);
        ClipPoint1.z = B.z + t * (A.z - B.z);
        ClipPoint1.w = B.w;
        float t2 = (-C.w - C.z) / (A.w - C.w + A.z - C.z);
        ClipPoint2.x = C.x + t2 * (A.x - C.x);
        ClipPoint2.y = C.y + t2 * (A.y - C.y);
        ClipPoint2.z = C.z + t2 * (A.z - C.z);
        ClipPoint2.w = C.w;

        //B and C clipped, repalce with A, ClipPoint1, ClipPoint2

        // cltr.tr1.points[0].gl_Position.x = A.x;
        // cltr.tr1.points[0].gl_Position.y = A.y;
        // cltr.tr1.points[0].gl_Position.z = A.z;
        // cltr.tr1.points[0].gl_Position.w = A.w;
        set_coords(cltr.tr1, 0, A.x, A.y, A.z, A.w);
        // for (int i = 0; i < maxAttributes; i++){
          // cltr.tr1.points[0].attributes[i] = A.attributes[i];
        // }
        set_attribs(A, cltr, 1, 0);
        interpolate_attribs(B, A, cltr, 1, 1, t);
        interpolate_attribs(C, A, cltr, 2, 2, t2);

        // std::cerr << cltr.tr1.points[0].attributes[3].u1 << std::endl;

        // cltr.tr1.points[1].gl_Position.x = ClipPoint1.x;
        // cltr.tr1.points[1].gl_Position.y = ClipPoint1.y;
        // cltr.tr1.points[1].gl_Position.z = ClipPoint1.z;
        // cltr.tr1.points[1].gl_Position.w = ClipPoint1.w;
        set_coords(cltr.tr1, 1, ClipPoint1.x, ClipPoint1.y, ClipPoint1.z, ClipPoint1.w);
        // for (int i = 0; i < maxAttributes; i++){
          // cltr.tr1.points[1].attributes[i] = B.attributes[i];
        // }

        // cltr.tr1.points[2].gl_Position.x = ClipPoint2.x;
        // cltr.tr1.points[2].gl_Position.y = ClipPoint2.y;
        // cltr.tr1.points[2].gl_Position.z = ClipPoint2.z;
        // cltr.tr1.points[2].gl_Position.w = ClipPoint2.w;
        set_coords(cltr.tr1, 2, ClipPoint2.x, ClipPoint2.y, ClipPoint2.z, ClipPoint2.w);
        // for (int i = 0; i < maxAttributes; i++){
          // cltr.tr1.points[2].attributes[i] = C.attributes[i];
        // }

        
      }
      
      else if (B.isOutside)
      {
        float t = (-A.w - A.z) / (B.w - A.w + B.z - A.z);
        ClipPoint1.x = A.x + t * (B.x - A.x);
        ClipPoint1.y = A.y + t * (B.y - A.y);
        ClipPoint1.z = A.z + t * (B.z - A.z);
        ClipPoint1.w = A.w;
        float t2 = (-C.w - C.z) / (B.w - C.w + B.z - C.z);
        ClipPoint2.x = C.x + t2 * (B.x - C.x);
        ClipPoint2.y = C.y + t2 * (B.y - C.y);
        ClipPoint2.z = C.z + t2 * (B.z - C.z);
        ClipPoint2.w = C.w;

        //A and C clipped, repalce with B, ClipPoint1, ClipPoint2
        // cltr.tr1.points[0].gl_Position.x = ClipPoint1.x;
        // cltr.tr1.points[0].gl_Position.y = ClipPoint1.y;
        // cltr.tr1.points[0].gl_Position.z = ClipPoint1.z;
        // cltr.tr1.points[0].gl_Position.w = ClipPoint1.w;
        set_coords(cltr.tr1, 0, ClipPoint1.x, ClipPoint1.y, ClipPoint1.z, ClipPoint1.w);
        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr1.points[0].attributes[i] = A.attributes[i];
        // }

        interpolate_attribs(A, B, cltr, 1, 0, t);
        set_attribs(B, cltr, 1, 1);
        interpolate_attribs(C, B, cltr, 1, 2, t2);

        // std::cerr << cltr.tr1.points[0].attributes[3].u1 << std::endl;

        // cltr.tr1.points[1].gl_Position.x = B.x;
        // cltr.tr1.points[1].gl_Position.y = B.y;
        // cltr.tr1.points[1].gl_Position.z = B.z;
        // cltr.tr1.points[1].gl_Position.w = B.w;
        set_coords(cltr.tr1, 1, B.x, B.y, B.z, B.w);
        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr1.points[1].attributes[i] = B.attributes[i];
        // }

        // cltr.tr1.points[2].gl_Position.x = ClipPoint2.x;
        // cltr.tr1.points[2].gl_Position.y = ClipPoint2.y;
        // cltr.tr1.points[2].gl_Position.z = ClipPoint2.z;
        // cltr.tr1.points[2].gl_Position.w = ClipPoint2.w;
        set_coords(cltr.tr1, 2, ClipPoint2.x, ClipPoint2.y, ClipPoint2.z, ClipPoint2.w);
        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr1.points[2].attributes[i] = C.attributes[i];
        // }


      }

      else if (C.isOutside)
      {
        float t = (-A.w - A.z) / (C.w - A.w + C.z - A.z);
        // std::cerr << "t: " << t << std::endl;
        ClipPoint1.x = A.x + t * (C.x - A.x);
        ClipPoint1.y = A.y + t * (C.y - A.y);
        ClipPoint1.z = A.z + t * (C.z - A.z);
        ClipPoint1.w = A.w;
        float t2 = (-B.w - B.z) / (C.w - B.w + C.z - B.z);
        // std::cerr << "t2: " << t2 << std::endl;
        ClipPoint2.x = B.x + t2 * (C.x - B.x);
        ClipPoint2.y = B.y + t2 * (C.y - B.y);
        ClipPoint2.z = B.z + t2 * (C.z - B.z);
        ClipPoint2.w = B.w;

        // float t = (-B.w - B.z) / (A.w - B.w + A.z - B.z);
        // ClipPoint1.x = B.x + t * (A.x - B.x);
        // ClipPoint1.y = B.y + t * (A.y - B.y);
        // ClipPoint1.z = B.z + t * (A.z - B.z);
        // ClipPoint1.w = B.w;
        // float t2 = (-C.w - C.z) / (A.w - C.w + A.z - C.z);
        // ClipPoint2.x = C.x + t2 * (A.x - C.x);
        // ClipPoint2.y = C.y + t2 * (A.y - C.y);
        // ClipPoint2.z = C.z + t2 * (A.z - C.z);
        // ClipPoint2.w = C.w;
        

        //A and B clipped, repalce with C, ClipPoint1, ClipPoint2

        // cltr.tr1.points[0].gl_Position.x = ClipPoint1.x;
        // cltr.tr1.points[0].gl_Position.y = ClipPoint1.y;
        // cltr.tr1.points[0].gl_Position.z = ClipPoint1.z;
        // cltr.tr1.points[0].gl_Position.w = ClipPoint1.w;
        set_coords(cltr.tr1, 0, ClipPoint1.x, ClipPoint1.y, ClipPoint1.z, ClipPoint1.w);
        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr1.points[0].attributes[i] = A.attributes[i];
        // }

        interpolate_attribs(A, C, cltr, 1, 0, t);
        interpolate_attribs(B, C, cltr, 1, 1, t2);
        set_attribs(C, cltr, 1, 2);

        // cltr.tr1.points[1].gl_Position.x = ClipPoint2.x;
        // cltr.tr1.points[1].gl_Position.y = ClipPoint2.y;
        // cltr.tr1.points[1].gl_Position.z = ClipPoint2.z;
        // cltr.tr1.points[1].gl_Position.w = ClipPoint2.w;
        set_coords(cltr.tr1, 1, ClipPoint2.x, ClipPoint2.y, ClipPoint2.z, ClipPoint2.w);
        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr1.points[1].attributes[i] = B.attributes[i];
        // }

        // cltr.tr1.points[2].gl_Position.x = C.x;
        // cltr.tr1.points[2].gl_Position.y = C.y;
        // cltr.tr1.points[2].gl_Position.z = C.z;
        // cltr.tr1.points[2].gl_Position.w = C.w;
        set_coords(cltr.tr1, 2, C.x, C.y, C.z, C.w);
        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr1.points[2].attributes[i] = C.attributes[i];
        // }

      }

      // cltr.tr1.


    }

    if (outside == 2 && inside == 1){
      if (A.isInside)
      {
        // float t = (-B.w - B.z) / (A.w - B.w + A.z - B.z);
        float t = (-A.w - A.z) / (C.w - A.w + C.z - A.z);
        ClipPoint1.x = A.x + t * (C.x - A.x);
        ClipPoint1.y = A.y + t * (C.y - A.y);
        ClipPoint1.z = A.z + t * (C.z - A.z);
        ClipPoint1.w = A.w;
        float t2 = (-A.w - A.z) / (B.w - A.w + B.z - A.z);
        ClipPoint2.x = A.x + t2 * (B.x - A.x);
        ClipPoint2.y = A.y + t2 * (B.y - A.y);
        ClipPoint2.z = A.z + t2 * (B.z - A.z);
        ClipPoint2.w = A.w;

        //B and C clipped, repalce with A, ClipPoint1, ClipPoint2
        interpolate_attribs(A, C, cltr, 1, 0, t);
        interpolate_attribs(A, B, cltr, 1, 1, t2);
        set_attribs(B, cltr, 1, 2);

        interpolate_attribs(A, C, cltr, 2, 0, t);
        set_attribs(B, cltr, 2, 1);
        set_attribs(C, cltr, 2, 2);
        // cltr.tr1.points[0].gl_Position.x = ClipPoint1.x;
        // cltr.tr1.points[0].gl_Position.y = ClipPoint1.y;
        // cltr.tr1.points[0].gl_Position.z = ClipPoint1.z;
        // cltr.tr1.points[0].gl_Position.w = ClipPoint1.w;

        // cltr.tr1.points[1].gl_Position.x = ClipPoint2.x;
        // cltr.tr1.points[1].gl_Position.y = ClipPoint2.y;
        // cltr.tr1.points[1].gl_Position.z = ClipPoint2.z;
        // cltr.tr1.points[1].gl_Position.w = ClipPoint2.w;

        // cltr.tr1.points[2].gl_Position.x = B.x;
        // cltr.tr1.points[2].gl_Position.y = B.y;
        // cltr.tr1.points[2].gl_Position.z = B.z;
        // cltr.tr1.points[2].gl_Position.w = B.w;

        set_coords(cltr.tr1, 0, ClipPoint1.x, ClipPoint1.y, ClipPoint1.z, ClipPoint1.w);
        set_coords(cltr.tr1, 1, ClipPoint2.x, ClipPoint2.y, ClipPoint2.z, ClipPoint2.w);
        set_coords(cltr.tr1, 2, B.x, B.y, B.z, B.w);
        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr1.points[0].attributes[i] = A.attributes[i];
        //   cltr.tr1.points[1].attributes[i] = A.attributes[i];
        //   cltr.tr1.points[2].attributes[i] = B.attributes[i];
        // }

        //------------

        // cltr.tr2.points[0].gl_Position.x = ClipPoint1.x;
        // cltr.tr2.points[0].gl_Position.y = ClipPoint1.y;
        // cltr.tr2.points[0].gl_Position.z = ClipPoint1.z;
        // cltr.tr2.points[0].gl_Position.w = ClipPoint1.w;

        set_coords(cltr.tr2, 0, ClipPoint1.x, ClipPoint1.y, ClipPoint1.z, ClipPoint1.w);
        set_coords(cltr.tr2, 1, B.x, B.y, B.z, B.w);
        set_coords(cltr.tr2, 2, C.x, C.y, C.z, C.w);
        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr2.points[0].attributes[i] = A.attributes[i];
        //   cltr.tr2.points[1].attributes[i] = B.attributes[i];
        //   cltr.tr2.points[2].attributes[i] = C.attributes[i];
        // }
        
      }

      if (B.isInside){
        // float t = (-C.w - C.z) / (B.w - C.w + B.z - C.z);
        float t = (-B.w - B.z) / (A.w - B.w + A.z - B.z);
        ClipPoint1.x = B.x + t * (A.x - B.x);
        ClipPoint1.y = B.y + t * (A.y - B.y);
        ClipPoint1.z = B.z + t * (A.z - B.z);
        ClipPoint1.w = B.w;
        float t2 = (-B.w - B.z) / (C.w - B.w + C.z - B.z);
        ClipPoint2.x = B.x + t2 * (C.x - B.x);
        ClipPoint2.y = B.y + t2 * (C.y - B.y);
        ClipPoint2.z = B.z + t2 * (C.z - B.z);
        ClipPoint2.w = B.w;

        //A and C clipped, repalce with B, ClipPoint1, ClipPoint2

        interpolate_attribs(B, A, cltr, 1, 0, t);
        interpolate_attribs(B, C, cltr, 1, 1, t2);
        set_attribs(C, cltr, 1, 2);

        interpolate_attribs(B, A, cltr, 2, 0, t);
        set_attribs(C, cltr, 2, 1);
        set_attribs(A, cltr, 2, 2);

        // cltr.tr1.points[0].gl_Position.x = ClipPoint1.x;
        // cltr.tr1.points[0].gl_Position.y = ClipPoint1.y;
        // cltr.tr1.points[0].gl_Position.z = ClipPoint1.z;
        // cltr.tr1.points[0].gl_Position.w = ClipPoint1.w;

        // cltr.tr1.points[1].gl_Position.x = ClipPoint2.x;
        // cltr.tr1.points[1].gl_Position.y = ClipPoint2.y;
        // cltr.tr1.points[1].gl_Position.z = ClipPoint2.z;
        // cltr.tr1.points[1].gl_Position.w = ClipPoint2.w;

        // cltr.tr1.points[2].gl_Position.x = C.x;
        // cltr.tr1.points[2].gl_Position.y = C.y;
        // cltr.tr1.points[2].gl_Position.z = C.z;
        // cltr.tr1.points[2].gl_Position.w = C.w;

        set_coords(cltr.tr1, 0, ClipPoint1.x, ClipPoint1.y, ClipPoint1.z, ClipPoint1.w);
        set_coords(cltr.tr1, 1, ClipPoint2.x, ClipPoint2.y, ClipPoint2.z, ClipPoint2.w);
        set_coords(cltr.tr1, 2, C.x, C.y, C.z, C.w);

        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr1.points[0].attributes[i] = B.attributes[i];
        //   cltr.tr1.points[1].attributes[i] = B.attributes[i];
        //   cltr.tr1.points[2].attributes[i] = C.attributes[i];
        // }

        //------------

        set_coords(cltr.tr2, 0, ClipPoint1.x, ClipPoint1.y, ClipPoint1.z, ClipPoint1.w);
        set_coords(cltr.tr2, 1, C.x, C.y, C.z, C.w);
        set_coords(cltr.tr2, 2, A.x, A.y, A.z, A.w);

        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr2.points[0].attributes[i] = B.attributes[i];
        //   cltr.tr2.points[1].attributes[i] = C.attributes[i];
        //   cltr.tr2.points[2].attributes[i] = A.attributes[i];
        // }
      }
      
      if (C.isInside){
        // float t = (-A.w - A.z) / (C.w - A.w + C.z - A.z);
        float t = (-C.w - C.z) / (B.w - C.w + B.z - C.z);
        ClipPoint1.x = C.x + t * (B.x - C.x);
        ClipPoint1.y = C.y + t * (B.y - C.y);
        ClipPoint1.z = C.z + t * (B.z - C.z);
        ClipPoint1.w = C.w;
        float t2 = (-C.w - C.z) / (A.w - C.w + A.z - C.z);
        ClipPoint2.x = C.x + t2 * (A.x - C.x);
        ClipPoint2.y = C.y + t2 * (A.y - C.y);
        ClipPoint2.z = C.z + t2 * (A.z - C.z);
        ClipPoint2.w = C.w;

        //B and A clipped, repalce with C, ClipPoint1, ClipPoint2

        interpolate_attribs(C, B, cltr, 1, 0, t);
        interpolate_attribs(C, A, cltr, 1, 1, t2);
        set_attribs(A, cltr, 1, 2);

        interpolate_attribs(C, B, cltr, 2, 0, t);
        set_attribs(A, cltr, 2, 1);
        set_attribs(B, cltr, 2, 2);

        // cltr.tr1.points[0].gl_Position.x = ClipPoint1.x;
        // cltr.tr1.points[0].gl_Position.y = ClipPoint1.y;
        // cltr.tr1.points[0].gl_Position.z = ClipPoint1.z;
        // cltr.tr1.points[0].gl_Position.w = ClipPoint1.w;

        // cltr.tr1.points[1].gl_Position.x = ClipPoint2.x;
        // cltr.tr1.points[1].gl_Position.y = ClipPoint2.y;
        // cltr.tr1.points[1].gl_Position.z = ClipPoint2.z;
        // cltr.tr1.points[1].gl_Position.w = ClipPoint2.w;

        // cltr.tr1.points[2].gl_Position.x = A.x;
        // cltr.tr1.points[2].gl_Position.y = A.y;
        // cltr.tr1.points[2].gl_Position.z = A.z;
        // cltr.tr1.points[2].gl_Position.w = A.w;

        set_coords(cltr.tr1, 0, ClipPoint1.x, ClipPoint1.y, ClipPoint1.z, ClipPoint1.w);
        set_coords(cltr.tr1, 1, ClipPoint2.x, ClipPoint2.y, ClipPoint2.z, ClipPoint2.w);
        set_coords(cltr.tr1, 2, A.x, A.y, A.z, A.w);

        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr1.points[0].attributes[i] = C.attributes[i];
        //   cltr.tr1.points[1].attributes[i] = C.attributes[i];
        //   cltr.tr1.points[2].attributes[i] = A.attributes[i];
        // }

        //------------

        set_coords(cltr.tr2, 0, ClipPoint1.x, ClipPoint1.y, ClipPoint1.z, ClipPoint1.w);
        set_coords(cltr.tr2, 1, A.x, A.y, A.z, A.w);
        set_coords(cltr.tr2, 2, B.x, B.y, B.z, B.w);

        // for (int i = 0; i < maxAttributes; i++){
        //   cltr.tr2.points[0].attributes[i] = C.attributes[i];
        //   cltr.tr2.points[1].attributes[i] = A.attributes[i];
        //   cltr.tr2.points[2].attributes[i] = B.attributes[i];
        // }
      }
    }
}

void reset_clip_triangles(ClippedTriangle &cltr){

  int uninitialized = 0x7F800000;

  for (int i = 0; i < 3; i++){
    cltr.tr1.points[i].gl_Position.x = uninitialized;
    cltr.tr1.points[i].gl_Position.y = uninitialized;
    cltr.tr1.points[i].gl_Position.z = uninitialized;
    cltr.tr1.points[i].gl_Position.w = uninitialized;
  }

  for (int i = 0; i < 3; i++){
    cltr.tr2.points[i].gl_Position.x = uninitialized;
    cltr.tr2.points[i].gl_Position.y = uninitialized;
    cltr.tr2.points[i].gl_Position.z = uninitialized;
    cltr.tr2.points[i].gl_Position.w = uninitialized;
  }

}

void drawTrianglesImpl(GPUMemory &mem, int i, DrawCommand &cmd, int drawID){
  // DrawCommand cmd = cb.commands[i].data.drawCommand;
  int nofVertices = cmd.nofVertices;

  int t = 0;
  bool rasterized = true;
  while ( t < nofVertices ) {
    Triangle triangle;
    
    // p("im dead?");

    loadTriangle(triangle, t, mem, cmd, drawID);

    // p("i live!");

    // std::cerr << triangle.points[0].attributes[3].u1 << std::endl;

    ClippedTriangle clippedTriangle;
    reset_clip_triangles(clippedTriangle);
    clipping(mem, clippedTriangle, triangle, 0, mem.framebuffer);

    if (clippedTriangle.tr1.points[0].gl_Position == glm::vec4(0x7F800000) ||
        clippedTriangle.tr1.points[1].gl_Position == glm::vec4(0x7F800000) ||
        clippedTriangle.tr1.points[2].gl_Position == glm::vec4(0x7F800000)){

        rasterized = false;
    }
    else
    {
      perspectiveDivision(clippedTriangle.tr1);

      // point 0:
      // 0
      // 0
      // -1.5
      // point 1:
      // -1.5
      // 0
      // -1.5
      // point 2:
      // -1
      // -1
      // 1
      // point 0:
      // 0
      // 0
      // -1.5
      // point 1:
      // -1
      // -1
      // 1
      // point 2:
      // 1
      // -1
      // 1

      // p("point 0:");
      // std::cerr << clippedTriangle.tr1.points[0].gl_Position.x << std::endl;
      // std::cerr << clippedTriangle.tr1.points[0].gl_Position.y << std::endl;
      // std::cerr << clippedTriangle.tr1.points[0].gl_Position.z << std::endl;
      // p("point 1:");
      // std::cerr << clippedTriangle.tr1.points[1].gl_Position.x << std::endl;
      // std::cerr << clippedTriangle.tr1.points[1].gl_Position.y << std::endl;
      // std::cerr << clippedTriangle.tr1.points[1].gl_Position.z << std::endl;
      // p("point 2:");
      // std::cerr << clippedTriangle.tr1.points[2].gl_Position.x << std::endl;
      // std::cerr << clippedTriangle.tr1.points[2].gl_Position.y << std::endl;
      // std::cerr << clippedTriangle.tr1.points[2].gl_Position.z << std::endl;
      // p("1 triangle!!!!!!!!");
      // std::cerr << "Perspective division done " << t << std::endl;
      viewportTrans(clippedTriangle.tr1, mem.framebuffer.width, mem.framebuffer.height);
      // std::cerr << "Viewport done " << t << std::endl;
      rasterize(mem.framebuffer, clippedTriangle.tr1, mem.programs[cmd.programID],
                mem, cmd.backfaceCulling);
      // std::cerr << "Rasterize done " << t << std::endl;

                // std::cerr << "aaa " << std::endl;
    }

    if (clippedTriangle.tr2.points[0].gl_Position == glm::vec4(0x7F800000) ||
        clippedTriangle.tr2.points[1].gl_Position == glm::vec4(0x7F800000) ||
        clippedTriangle.tr2.points[2].gl_Position == glm::vec4(0x7F800000)){

        rasterized = false;
    }
    else
    {
      perspectiveDivision(clippedTriangle.tr2);

      // p("point 0:");
      // std::cerr << clippedTriangle.tr2.points[0].gl_Position.x << std::endl;
      // std::cerr << clippedTriangle.tr2.points[0].gl_Position.y << std::endl;
      // std::cerr << clippedTriangle.tr2.points[0].gl_Position.z << std::endl;
      // p("point 1:");
      // std::cerr << clippedTriangle.tr2.points[1].gl_Position.x << std::endl;
      // std::cerr << clippedTriangle.tr2.points[1].gl_Position.y << std::endl;
      // std::cerr << clippedTriangle.tr2.points[1].gl_Position.z << std::endl;
      // p("point 2:");
      // std::cerr << clippedTriangle.tr2.points[2].gl_Position.x << std::endl;
      // std::cerr << clippedTriangle.tr2.points[2].gl_Position.y << std::endl;
      // std::cerr << clippedTriangle.tr2.points[2].gl_Position.z << std::endl;

      // p("2 triangle!!!!!!!!");
      // std::cerr << "Perspective division done " << t << std::endl;
      viewportTrans(clippedTriangle.tr2, mem.framebuffer.width, mem.framebuffer.height);
      // std::cerr << "Viewport done " << t << std::endl;
      rasterize(mem.framebuffer, clippedTriangle.tr2, mem.programs[cmd.programID],
                mem, cmd.backfaceCulling);
      // std::cerr << "Rasterize done " << t << std::endl;

                // std::cerr << "bbbb " << std::endl;

    }

    // if (rasterized == false)
      // p("triangle skipped");

    // if (rasterized == true)
      // p("triangle ready");

    t += 3;

    // perspectiveDivision(triangle);
    // viewportTrans(triangle, mem.framebuffer.width, mem.framebuffer.height);
    // // std::cout << "triangle.points[0].gl_Position: " << triangle.points[0].gl_Position.x << " " << triangle.points[0].gl_Position.y << " " << triangle.points[0].gl_Position.z << " " << triangle.points[0].gl_Position.w << std::endl;
    // // viewportTransformation(triangle, ctx.frame.width, ctx.frame.height);
    // rasterize(mem.framebuffer, triangle, mem.programs[cmd.programID], mem, cmd.backfaceCulling);
    // std::cerr << "test1" << std::endl;
  }

  // p("end?????????????????????????");
}

// ! [gpu_execute]
void gpu_execute(GPUMemory &mem,CommandBuffer &cb){

  int drawID = -1;
  for(uint32_t i = 0; i < cb.nofCommands; i++){
    CommandType type = cb.commands[i].type;
    CommandData data = cb.commands[i].data;
    if(type == CommandType::CLEAR){
      clear(mem,data.clearCommand);
    }

    if(type == CommandType::DRAW){
      drawID++;
      // draw(mem, i, cb, drawID);
      drawTrianglesImpl(mem, i, data.drawCommand, drawID);
    }
  }

}
//! [gpu_execute]

/**
 * @brief This function reads color from texture.
 *
 * @param texture texture
 * @param uv uv coordinates
 *
 * @return color 4 floats
 */
glm::vec4 read_texture(Texture const&texture,glm::vec2 uv){

  if (!texture.data) return glm::vec4(0.f);

  auto uv1 = glm::fract(uv);
  auto uv2 = uv1 * glm::vec2(texture.width-1, texture.height-1) + 0.5f;
  auto pix = glm::uvec2(uv2);

  //auto t   = glm::fract(uv2);
  glm::vec4 color = glm::vec4(0.f, 0.f, 0.f, 1.f);

  for(uint32_t c = 0; c < texture.channels; ++c)
    color[c] = texture.data[(pix.y*texture.width + pix.x) * texture.channels + c] / 255.f;

  return color;
}

