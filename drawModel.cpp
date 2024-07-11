/*!
 * @file
 * @brief This file contains functions for model rendering
 *
 * @author Tomáš Milet, imilet@fit.vutbr.cz
 */
#include <student/drawModel.hpp>
#include <student/gpu.hpp>

#include <iostream>
#include <string>

///\endcond

void drawModel_vertexShader(OutVertex&outVertex,InVertex const&inVertex,ShaderInterface const&si);
void drawModel_fragmentShader(OutFragment&outFragment, InFragment const&inFragment, ShaderInterface const&si);


void pf(std::string s){
  std::cout << s << std::endl;
}

void clear(CommandBuffer &cb){
  cb.commands[0].type = CommandType::CLEAR;
  cb.commands[0].data.clearCommand.color = {0.1f, 0.15f, 0.1f, 1.0f};
  cb.commands[0].data.clearCommand.depth = 1e+11;
  cb.commands[0].data.clearCommand.clearColor = true;
  cb.commands[0].data.clearCommand.clearDepth = true;
  cb.nofCommands++;
}

void addDrawCommand(CommandBuffer &cb, bool BFCulling,
 IndexType indType, size_t indOffset, int32_t vaoID,
 size_t nofVertices,
 VertexAttrib pos, VertexAttrib norm, VertexAttrib tex){
  cb.commands[cb.nofCommands].type = CommandType::DRAW;
  cb.commands[cb.nofCommands].data.drawCommand.backfaceCulling = BFCulling;
  cb.commands[cb.nofCommands].data.drawCommand.vao.indexType = indType;
  cb.commands[cb.nofCommands].data.drawCommand.vao.indexOffset = indOffset;
  cb.commands[cb.nofCommands].data.drawCommand.vao.indexBufferID = vaoID;
  cb.commands[cb.nofCommands].data.drawCommand.nofVertices = nofVertices;
  cb.commands[cb.nofCommands].data.drawCommand.programID = 0;

  VertexAttrib *attribs[3] = {&pos, &norm, &tex};

  if (attribs[0] == nullptr){
    pf("ERROR: addDrawCommand: attribs[0] == nullptr");
    return;
  }

  if (attribs[1] == nullptr){
    pf("ERROR: addDrawCommand: attribs[1] == nullptr");
    return;
  }

  if (attribs[2] == nullptr){
    pf("ERROR: addDrawCommand: attribs[2] == nullptr");
    return;
  }

  for (int i = 0; i < maxAttributes-1; i++){

    if (attribs[i] == nullptr){
      cb.commands[cb.nofCommands].data.drawCommand.vao.vertexAttrib[i].bufferID = -1;
      continue;
    }

    cb.commands[cb.nofCommands].data.drawCommand.vao.vertexAttrib[i].bufferID = attribs[i]->bufferID;
    cb.commands[cb.nofCommands].data.drawCommand.vao.vertexAttrib[i].type = attribs[i]->type;
    cb.commands[cb.nofCommands].data.drawCommand.vao.vertexAttrib[i].stride = attribs[i]->stride;
    cb.commands[cb.nofCommands].data.drawCommand.vao.vertexAttrib[i].offset = attribs[i]->offset;
  }

  cb.nofCommands++;
}

void writeToMemory(GPUMemory&mem, Model const&model){
  mem.programs[0].vertexShader = drawModel_vertexShader;
  mem.programs[0].fragmentShader = drawModel_fragmentShader;

  mem.programs[0].vs2fs[0] = AttributeType::VEC3;
  mem.programs[0].vs2fs[1] = AttributeType::VEC3;
  mem.programs[0].vs2fs[2] = AttributeType::VEC2;
  mem.programs[0].vs2fs[3] = AttributeType::UINT;

  int modeltex_size = model.textures.size();

  for (int i = 0; i < modeltex_size; i++){
    mem.textures[i].data = model.textures[i].data;
    mem.textures[i].width = model.textures[i].width;
    mem.textures[i].height = model.textures[i].height;
    mem.textures[i].channels = model.textures[i].channels;
  }

  int modelbuf_size = model.buffers.size();

  for (int i = 0; i < modelbuf_size; i++){
    mem.buffers[i].data = model.buffers[i].data;
    mem.buffers[i].size = model.buffers[i].size;
  }

}

void SetMeshes(GPUMemory&mem, Model const&model, uint32_t MeshID,
               uint32_t DrawID, glm::mat4 Mat, glm::mat4 iMat){
  // ...
  // mem.programs[0].vertexShader = drawModel_vertexShader;
  // mem.programs[0].fragmentShader = drawModel_fragmentShader;


    mem.uniforms[10+DrawID*5+0].m4 = Mat;

    mem.uniforms[10+DrawID*5+1].m4 = iMat;

    glm::vec4 DifCOl = model.meshes[MeshID].diffuseColor;
    mem.uniforms[10+DrawID*5+2].v4 = DifCOl;

    uint32_t DifTex = model.meshes[MeshID].diffuseTexture;
    mem.uniforms[10+DrawID*5+3].i1 = DifTex;

    bool DbSided = model.meshes[MeshID].doubleSided;
    mem.uniforms[10+DrawID*5+4].v1 = DbSided;

}

void prepareNode(GPUMemory&mem, CommandBuffer&cb, Node const&node,
                Model const&model, glm::mat4 prubeznaMatice,
                glm::mat4 IpMat, 
                int32_t meshID, uint32_t &DrawID){
  
  if(node.mesh >= 0){
 
    Mesh const&mesh = model.meshes[node.mesh];
    
    //tvorba kreliciho prikazu
    // mesh.doubleSided ...;
    // mesh.position ...;
    // mesh.normal ...;
    // mesh.texCoord ...;
    // mesh.indexBufferID ...;
    // mesh.indexOffset ...;
    // mesh.indexType ...
    // mesh.nofVertices ...

    bool bfCulling = mesh.doubleSided;
    IndexType indType = mesh.indexType;
    size_t indOffset = mesh.indexOffset;
    int32_t vaoID = mesh.indexBufferID;
    uint32_t NumOVert = mesh.nofIndices;

    VertexAttrib pos = mesh.position;
    VertexAttrib norm = mesh.normal;
    VertexAttrib tex = mesh.texCoord;
 
    addDrawCommand(cb, !bfCulling, indType, indOffset, vaoID, NumOVert,
                  pos, norm, tex);

    SetMeshes(mem, model, meshID, DrawID, prubeznaMatice, IpMat);
    DrawID++;

  }
 
  for(size_t i = 0; i < node.children.size(); ++i){
    int32_t meshI = node.children[i].mesh;
    glm::mat4 mat = prubeznaMatice * node.children[i].modelMatrix;

    glm::mat4 imat;
    imat = glm::transpose(glm::inverse(mat));
    
    prepareNode(mem, cb, node.children[i], model, mat, imat, meshI, DrawID);
  }
    //  rekurze
}



/**
 * @brief This function prepares model into memory and creates command buffer
 *
 * @param mem gpu memory
 * @param commandBuffer command buffer
 * @param model model structure
 */
//! [drawModel]
void prepareModel(GPUMemory&mem,CommandBuffer&commandBuffer,Model const&model){
  // (void)mem;
  // (void)commandBuffer;
  // (void)model;
  /// \todo Tato funkce připraví command buffer pro model a nastaví správně pamět grafické karty.<br>
  /// Vaším úkolem je správně projít model a vložit vykreslovací příkazy do commandBufferu.
  /// Zároveň musíte vložit do paměti textury, buffery a uniformní proměnné, které buffer command buffer využívat.
  /// Bližší informace jsou uvedeny na hlavní stránce dokumentace a v testech.

  CommandBuffer& cb = commandBuffer;

  clear(cb);

  uint32_t DrawID = 0;
  // model.roots[0]

  for(size_t i = 0; i < model.roots.size(); ++i){
      uint32_t beginID = model.roots[i].mesh;
      glm::mat4 mat = model.roots[i].modelMatrix;
      glm::mat4 Imat = glm::inverse(mat);
      prepareNode(mem, cb, model.roots[i], model, mat, Imat, beginID, DrawID);
  }
  writeToMemory(mem, model);

}
//! [drawModel]
void printMatrix(const glm::mat4& matrix) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      std::cout << matrix[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}
/**
 * @brief This function represents vertex shader of texture rendering method.
 *
 * @param outVertex output vertex
 * @param inVertex input vertex
 * @param si shader interface
 */
//! [drawModel_vs]
void drawModel_vertexShader(OutVertex& outVertex, InVertex const& inVertex, ShaderInterface const& si) {
  auto const pos = glm::vec4(inVertex.attributes[0].v3, 1.0f);
  auto const& nor = inVertex.attributes[1].v3;

  if (si.uniforms == nullptr)
    return;

  auto const&viewMatrix       = si.uniforms[0].m4;
  auto const&projectionMatrix = si.uniforms[1].m4;

    auto mvp = projectionMatrix*viewMatrix;

    auto const& modelMatrix = si.uniforms[10 + inVertex.gl_DrawID * 5 + 0].m4;
  auto const& invModelMatrix = glm::transpose(glm::inverse(modelMatrix));

  auto worldPosition = modelMatrix * pos;
  auto clipPosition = mvp * worldPosition;

  auto worldNormal = glm::vec3(invModelMatrix * glm::vec4(nor, 0.0f));

  auto texture = inVertex.attributes[2].v2;
  auto drawid = inVertex.gl_DrawID;

  // pf("hi4");
  outVertex.gl_Position = clipPosition;
  outVertex.attributes[0].v3 = worldPosition;
  outVertex.attributes[1].v3 = worldNormal;
  outVertex.attributes[2].v2 = texture;
  outVertex.attributes[3].u1 = drawid;
  // pf("hi5");
}
//! [drawModel_vs]

/**
 * @brief This functionrepresents fragment shader of texture rendering method.
 *
 * @param outFragment output fragment
 * @param inFragment input fragment
 * @param si shader interface
 */
//! [drawModel_fs]
void drawModel_fragmentShader(OutFragment&outFragment, InFragment const&inFragment, ShaderInterface const&si){
  // (void)outFragment;
  // (void)inFragment;
  // (void)si;
  /// \todo Tato funkce reprezentujte fragment shader.<br>
  /// Vaším úkolem je správně obarvit fragmenty a osvětlit je pomocí lambertova osvětlovacího modelu.
  /// Bližší informace jsou uvedeny na hlavní stránce dokumentace.
  // auto const& light          = si.uniforms[2].v3;
  // auto const& cameraPosition = si.uniforms[3].v3;
  // auto const& vpos           = inFragment.attributes[0].v3;
  // auto const& vnor           = inFragment.attributes[1].v3;
  // auto vvnor = glm::normalize(vnor);
  // outFragment.gl_FragColor = glm::vec4(color,1.f);

  auto position = inFragment.attributes[0].v3;
  auto normal = inFragment.attributes[1].v3;
  auto texture_xy = inFragment.attributes[2].v2;
  auto drawid = inFragment.attributes[3].u1;

  // std::cout << "drawid: " << drawid << std::endl;

  auto light_pos = si.uniforms[1].v3;
  auto cam_pos = si.uniforms[2].v3;

  auto diffuse_color = si.uniforms[10 + drawid * 5 + 2].v4;
  auto texid = si.uniforms[10 + drawid * 5 + 3].i1;

  glm::vec2 texcoord = glm::vec2(0.f);
  Texture tex;

  glm::vec4 tex_color = glm::vec4(-1.f);

  if (texid != -1){
    tex = si.textures[texid];
    texcoord = glm::vec2(texture_xy.x, texture_xy.y);
    tex_color = read_texture(tex, texcoord);
  }
  auto doublesided = si.uniforms[10 + drawid * 5 + 4].v1;
  auto CamToSurface = glm::normalize(position - cam_pos);



  if (doublesided == 1.f && glm::dot(CamToSurface, glm::normalize(normal)) >= 0.f){
    normal = -normal;
  }

  glm::vec4 diffuse = diffuse_color;



  auto diffuse_factor = glm::clamp(glm::dot(glm::normalize(light_pos - position), glm::normalize(normal)), 0.f, 1.f);
  // glm::clamp(glm::dot(L,N),0.f,1.f)
  auto ambient_factor = 0.2f;

  if (tex_color != glm::vec4(-1.f)){
    diffuse = tex_color;
  }

  auto ambient_component = diffuse * ambient_factor;
  auto diffuse_component = diffuse * diffuse_factor;

  auto final_color = glm::vec4(ambient_component + diffuse_component);
  final_color.w = diffuse.w;

  outFragment.gl_FragColor = final_color;
}
//! [drawModel_fs]

